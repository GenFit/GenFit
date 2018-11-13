/* Copyright 2008-2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "KalmanFitterInfo.h"

#include <cassert>
#include <TBuffer.h>

#include "IO.h"
#include "Exception.h"
#include "FitStatus.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

//#define DEBUG


namespace genfit {

KalmanFitterInfo::KalmanFitterInfo() :
  AbsFitterInfo(), fixWeights_(false)
{
  ;
}

KalmanFitterInfo::KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)  :
  AbsFitterInfo(trackPoint, rep), fixWeights_(false)
{
  ;
}

KalmanFitterInfo::~KalmanFitterInfo() {
  deleteMeasurementInfo();
}


KalmanFitterInfo* KalmanFitterInfo::clone() const {
  KalmanFitterInfo* retVal = new KalmanFitterInfo(this->getTrackPoint(), this->getRep());
  if (hasReferenceState())
    retVal->setReferenceState(new ReferenceStateOnPlane(*getReferenceState()));
  if (hasForwardPrediction())
    retVal->setForwardPrediction(new MeasuredStateOnPlane(*getForwardPrediction()));
  if (hasForwardUpdate())
    retVal->setForwardUpdate(new KalmanFittedStateOnPlane(*getForwardUpdate()));
  if (hasBackwardPrediction())
    retVal->setBackwardPrediction(new MeasuredStateOnPlane(*getBackwardPrediction()));
  if (hasBackwardUpdate())
    retVal->setBackwardUpdate(new KalmanFittedStateOnPlane(*getBackwardUpdate()));

  retVal->measurementsOnPlane_.reserve(getNumMeasurements());
  for (std::vector<MeasurementOnPlane*>::const_iterator it = this->measurementsOnPlane_.begin(); it != this->measurementsOnPlane_.end(); ++it) {
    retVal->addMeasurementOnPlane(new MeasurementOnPlane(**it));
  }

  retVal->fixWeights_ = this->fixWeights_;

  return retVal;
}


MeasurementOnPlane KalmanFitterInfo::getAvgWeightedMeasurementOnPlane(bool ignoreWeights) const {

  MeasurementOnPlane retVal(*(measurementsOnPlane_[0]));

  if(measurementsOnPlane_.size() == 1) {
    if (ignoreWeights) {
      retVal.setWeight(1.);
    }
    else {
      double weight = (measurementsOnPlane_[0])->getWeight();
      if (weight != 1.) {
        retVal.getCov() *= 1. / weight;
      }
      retVal.setWeight(weight);
    }
    return retVal;
  }

  // more than one hit

  double sumOfWeights(0), weight(0);

  retVal.getState().Zero();
  retVal.getCov().Zero();

  TMatrixDSym covInv;

  for(unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {

    if (i>0) {
      // make sure we have compatible measurement types
      // TODO: replace with Exceptions!
      assert(measurementsOnPlane_[i]->getPlane() == measurementsOnPlane_[0]->getPlane());
      assert(*(measurementsOnPlane_[i]->getHMatrix()) == *(measurementsOnPlane_[0]->getHMatrix()));
    }

    tools::invertMatrix(measurementsOnPlane_[i]->getCov(), covInv); // invert cov
    if (ignoreWeights) {
      sumOfWeights += 1.;
    }
    else {
      weight = measurementsOnPlane_[i]->getWeight();
      sumOfWeights += weight;
      covInv *= weight; // weigh cov
    }
    retVal.getCov() += covInv; // covInv is already inverted and weighted

    retVal.getState() += covInv * measurementsOnPlane_[i]->getState();
  }

  // invert Cov
  tools::invertMatrix(retVal.getCov());

  retVal.getState() *= retVal.getCov();

  retVal.setWeight(sumOfWeights);

  return retVal;
}


MeasurementOnPlane* KalmanFitterInfo::getClosestMeasurementOnPlane(const StateOnPlane* sop) const {
  if(measurementsOnPlane_.size() == 0)
    return nullptr;

  if(measurementsOnPlane_.size() == 1)
    return getMeasurementOnPlane(0);

  double normMin(9.99E99);
  unsigned int iMin(0);
  const AbsHMatrix* H = measurementsOnPlane_[0]->getHMatrix();
  for (unsigned int i=0; i<getNumMeasurements(); ++i) {
    if (*(measurementsOnPlane_[i]->getHMatrix()) != *H){
      Exception e("KalmanFitterInfo::getClosestMeasurementOnPlane: Cannot compare measurements with different H-Matrices. Maybe you have to select a different multipleMeasurementHandling.", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }

    TVectorD res = measurementsOnPlane_[i]->getState() - H->Hv(sop->getState());
    double norm = sqrt(res.Norm2Sqr());
    if (norm < normMin) {
      normMin = norm;
      iMin = i;
    }
  }

  return getMeasurementOnPlane(iMin);
}


std::vector<double> KalmanFitterInfo::getWeights() const {
  std::vector<double> retVal(getNumMeasurements());

  for (unsigned int i=0; i<getNumMeasurements(); ++i) {
    retVal[i] = getMeasurementOnPlane(i)->getWeight();
  }

  return retVal;
}


const MeasuredStateOnPlane& KalmanFitterInfo::getFittedState(bool biased) const {

  // check if we can use cached states
  if (biased && fittedStateBiased_)
    return *fittedStateBiased_;
  if (!biased && fittedStateUnbiased_)
    return *fittedStateUnbiased_;


  const TrackPoint* tp = this->getTrackPoint();
  const Track* tr = tp->getTrack();
  const AbsTrackRep* rep =  this->getRep();

  bool first(false), last(false);
  PruneFlags& flag = tr->getFitStatus(rep)->getPruneFlags();
  // if Track is pruned so that only one TrackPoint remains, see if it was the first or last one
  #ifdef DEBUG
  if (flag.isPruned()) {
    debugOut << "KalmanFitterInfo::getFittedState - Track is pruned and has " << tr->getNumPoints() << " TrackPoints \n";
    flag.Print();
  }
  #endif
  if (flag.isPruned() && tr->getNumPoints() == 1) {
    if (flag.hasFlags("F")) {
      first = true;
      #ifdef DEBUG
      debugOut << "KalmanFitterInfo::getFittedState - has flag F \n";
      #endif
    }
    else if (flag.hasFlags("L")) {
      last = true;
      #ifdef DEBUG
      debugOut << "KalmanFitterInfo::getFittedState - has flag L \n";
      #endif
    }
  }
  else { // otherwise check against TrackPoint order
    first = tr->getPointWithFitterInfo(0, rep) == tp;
    last = tr->getPointWithFitterInfo(-1, rep) == tp;
  }

  #ifdef DEBUG
  debugOut << "KalmanFitterInfo::getFittedState first " << first << ", last " << last << "\n";
  debugOut << "KalmanFitterInfo::getFittedState forwardPrediction_ " << forwardPrediction_.get() << ", forwardUpdate_ " << forwardUpdate_.get() << "\n";
  debugOut << "KalmanFitterInfo::getFittedState backwardPrediction_ " << backwardPrediction_.get() << ", backwardUpdate_ " << backwardUpdate_.get() << "\n";
  #endif

  /*
  if both needed forward prediction/update and backward prediction are available,
  use them to calculate smoothed state.
  Otherwise, if one is missing (i.e. has been pruned) and we are at first or last hit,
  use only backward or forward prediction (unbiased) of update (biased).
  */

  if (biased) {
    // last measurement and no backward prediction -> use forward update
    if (last && !backwardPrediction_) {
      if(!forwardUpdate_) {
        Exception e("KalmanFitterInfo::getFittedState: Needed updates/predictions not available in this FitterInfo.", __LINE__,__FILE__);
        e.setFatal();
        throw e;
      }
      #ifdef DEBUG
      debugOut << "KalmanFitterInfo::getFittedState - biased at last measurement = forwardUpdate_ \n";
      #endif
      return *forwardUpdate_;
    }

    // first measurement and no forward update -> use backward update
    if (first && (!forwardUpdate_ || (backwardUpdate_ && !forwardPrediction_) ) ) {
      if(!backwardUpdate_.get()) {
        Exception e("KalmanFitterInfo::getFittedState: Needed updates/predictions not available in this FitterInfo.", __LINE__,__FILE__);
        e.setFatal();
        throw e;
      }
      #ifdef DEBUG
      debugOut << "KalmanFitterInfo::getFittedState - biased at first measurement = backwardUpdate_ \n";
      //backwardUpdate_->Print();
      #endif
      return *backwardUpdate_;
    }

    if(!forwardUpdate_ || !backwardPrediction_) {
      Exception e("KalmanFitterInfo::getFittedState: Needed updates/predictions not available in this FitterInfo.", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    #ifdef DEBUG
    debugOut << "KalmanFitterInfo::getFittedState - biased = mean(forwardUpdate_, backwardPrediction_) \n";
    #endif
    fittedStateBiased_.reset(new MeasuredStateOnPlane(calcAverageState(*forwardUpdate_, *backwardPrediction_)));
    return *fittedStateBiased_;
  }

  // unbiased

  // last measurement and no backward prediction -> use forward prediction
  if (last && !backwardPrediction_) {
    if(!forwardPrediction_) {
      Exception e("KalmanFitterInfo::getFittedState: Needed updates/predictions not available in this FitterInfo.", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    #ifdef DEBUG
    debugOut << "KalmanFitterInfo::getFittedState - unbiased at last measurement = forwardPrediction_ \n";
    #endif
    return *forwardPrediction_;
  }

  // first measurement and no forward prediction -> use backward prediction
  if (first && !forwardPrediction_) {
    if(!backwardPrediction_) {
      Exception e("KalmanFitterInfo::getFittedState: Needed updates/predictions not available in this FitterInfo.", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    #ifdef DEBUG
    debugOut << "KalmanFitterInfo::getFittedState - unbiased at first measurement = backwardPrediction_ \n";
    #endif
    return *backwardPrediction_;
  }

  if(!forwardPrediction_ || !backwardPrediction_) {
    Exception e("KalmanFitterInfo::getFittedState: Needed updates/predictions not available in this FitterInfo.", __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  #ifdef DEBUG
  debugOut << "KalmanFitterInfo::getFittedState - unbiased = mean(forwardPrediction_, backwardPrediction_) \n";
  #endif
  fittedStateUnbiased_.reset(new MeasuredStateOnPlane(calcAverageState(*forwardPrediction_, *backwardPrediction_)));
  return *fittedStateUnbiased_;
}


MeasurementOnPlane KalmanFitterInfo::getResidual(unsigned int iMeasurement, bool biased, bool onlyMeasurementErrors) const {

  const MeasuredStateOnPlane& smoothedState = getFittedState(biased);
  const MeasurementOnPlane* measurement = measurementsOnPlane_.at(iMeasurement);
  const SharedPlanePtr& plane = measurement->getPlane();

  // check equality of planes and reps
  if(*(smoothedState.getPlane()) != *plane) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }
  if(smoothedState.getRep() != measurement->getRep()) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined wrt the same TrackRep.", __LINE__,__FILE__);
    throw e;
  }

  const AbsHMatrix* H = measurement->getHMatrix();

  //TODO: shouldn't the definition be (smoothed - measured) ?
  // res = -(H*smoothedState - measuredState)
  TVectorD res(H->Hv(smoothedState.getState()));
  res -= measurement->getState();
  res *= -1;

  if (onlyMeasurementErrors) {
    return MeasurementOnPlane(res, measurement->getCov(), plane, smoothedState.getRep(), H->clone(), measurement->getWeight());
  }
    
  TMatrixDSym cov(smoothedState.getCov());
  H->HMHt(cov);
  cov += measurement->getCov();

  return MeasurementOnPlane(res, cov, plane, smoothedState.getRep(), H->clone(), measurement->getWeight());
}


double KalmanFitterInfo::getSmoothedChi2(unsigned int iMeasurement) const {
  const MeasurementOnPlane& res = getResidual(iMeasurement, true, false);

  TMatrixDSym Rinv;
  tools::invertMatrix(res.getCov(), Rinv);
  return Rinv.Similarity(res.getState());
}


void KalmanFitterInfo::setReferenceState(ReferenceStateOnPlane* referenceState) {
  referenceState_.reset(referenceState);
  if (referenceState_)
    setPlane(referenceState_->getPlane());

  // if plane has changed, delete outdated info
  /*if (forwardPrediction_ && forwardPrediction_->getPlane() != getPlane())
    setForwardPrediction(0);

  if (forwardUpdate_ && forwardUpdate_->getPlane() != getPlane())
    setForwardUpdate(0);

  if (backwardPrediction_ && backwardPrediction_->getPlane() != getPlane())
    setBackwardPrediction(0);

  if (backwardUpdate_ && backwardUpdate_->getPlane() != getPlane())
    setBackwardUpdate(0);

  if (measurementsOnPlane_.size() > 0 && measurementsOnPlane_[0]->getPlane() != getPlane())
    deleteMeasurementInfo();
  */
}


void KalmanFitterInfo::setForwardPrediction(MeasuredStateOnPlane* forwardPrediction) {
  forwardPrediction_.reset(forwardPrediction);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (forwardPrediction_)
    setPlane(forwardPrediction_->getPlane());
}

void KalmanFitterInfo::setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction) {
  backwardPrediction_.reset(backwardPrediction);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (backwardPrediction_)
    setPlane(backwardPrediction_->getPlane());
}

void KalmanFitterInfo::setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate) {
  forwardUpdate_.reset(forwardUpdate);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (forwardUpdate_)
    setPlane(forwardUpdate_->getPlane());
}

void KalmanFitterInfo::setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate) {
  backwardUpdate_.reset(backwardUpdate);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
  if (backwardUpdate_)
    setPlane(backwardUpdate_->getPlane());
}


void KalmanFitterInfo::setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane) {
  deleteMeasurementInfo();

  for (std::vector<MeasurementOnPlane*>::const_iterator m = measurementsOnPlane.begin(), mend = measurementsOnPlane.end(); m < mend; ++m) {
    addMeasurementOnPlane(*m);
  }
}


void KalmanFitterInfo::addMeasurementOnPlane(MeasurementOnPlane* measurementOnPlane) {
  if (measurementsOnPlane_.size() == 0)
    setPlane(measurementOnPlane->getPlane());

  measurementsOnPlane_.push_back(measurementOnPlane);
}

void KalmanFitterInfo::addMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane) {
  for (std::vector<MeasurementOnPlane*>::const_iterator m = measurementsOnPlane.begin(), mend = measurementsOnPlane.end(); m < mend; ++m) {
    addMeasurementOnPlane(*m);
  }
}


void KalmanFitterInfo::setRep(const AbsTrackRep* rep) {
  rep_ = rep;

  if (referenceState_)
    referenceState_->setRep(rep);

  if (forwardPrediction_)
    forwardPrediction_->setRep(rep);

  if (forwardUpdate_)
    forwardUpdate_->setRep(rep);

  if (backwardPrediction_)
    backwardPrediction_->setRep(rep);

  if (backwardUpdate_)
    backwardUpdate_->setRep(rep);

  for (std::vector<MeasurementOnPlane*>::iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    (*it)->setRep(rep);
  }
}


void KalmanFitterInfo::setWeights(const std::vector<double>& weights) {

  if (weights.size() != getNumMeasurements()) {
    Exception e("KalmanFitterInfo::setWeights: weights do not have the same size as measurementsOnPlane", __LINE__,__FILE__);
    throw e;
  }

  if (fixWeights_) {
    errorOut << "KalmanFitterInfo::setWeights - WARNING: setting weights even though weights are fixed!" << std::endl;
  }

  for (unsigned int i=0; i<getNumMeasurements(); ++i) {
    getMeasurementOnPlane(i)->setWeight(weights[i]);
  }
}


void KalmanFitterInfo::deleteForwardInfo() {
  setForwardPrediction(nullptr);
  setForwardUpdate(nullptr);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
}

void KalmanFitterInfo::deleteBackwardInfo() {
  setBackwardPrediction(nullptr);
  setBackwardUpdate(nullptr);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
}

void KalmanFitterInfo::deletePredictions() {
  setForwardPrediction(nullptr);
  setBackwardPrediction(nullptr);
  fittedStateUnbiased_.reset();
  fittedStateBiased_.reset();
}

void KalmanFitterInfo::deleteMeasurementInfo() {
  // FIXME: need smart pointers / smart containers here
  for (size_t i = 0; i < measurementsOnPlane_.size(); ++i)
    delete measurementsOnPlane_[i];

  measurementsOnPlane_.clear();
}


void KalmanFitterInfo::Print(const Option_t*) const {
  printOut << "genfit::KalmanFitterInfo. Belongs to TrackPoint " << trackPoint_ << "; TrackRep " <<  rep_  << "\n";

  if (fixWeights_)
    printOut << "Weights are fixed.\n";

  for (unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {
    printOut << "MeasurementOnPlane Nr " << i <<": "; measurementsOnPlane_[i]->Print();
  }

  if (referenceState_) {
    printOut << "Reference state: "; referenceState_->Print();
  }
  if (forwardPrediction_) {
    printOut << "Forward prediction_: "; forwardPrediction_->Print();
  }
  if (forwardUpdate_) {
    printOut << "Forward update: "; forwardUpdate_->Print();
  }
  if (backwardPrediction_) {
    printOut << "Backward prediction_: "; backwardPrediction_->Print();
  }
  if (backwardUpdate_) {
    printOut << "Backward update: "; backwardUpdate_->Print();
  }

}


bool KalmanFitterInfo::checkConsistency(const genfit::PruneFlags* flags) const {

  bool retVal(true);

  // check if in a TrackPoint
  if (!trackPoint_) {
    errorOut << "KalmanFitterInfo::checkConsistency(): trackPoint_ is nullptr" << std::endl;
    retVal = false;
  }

  // check if there is a reference state
  /*if (!referenceState_) {
    errorOut << "KalmanFitterInfo::checkConsistency(): referenceState_ is nullptr" << std::endl;
    retVal = false;
  }*/

  SharedPlanePtr plane = getPlane();

  if (plane.get() == nullptr) {
    if (!(referenceState_ || forwardPrediction_ || forwardUpdate_ || backwardPrediction_ || backwardUpdate_ || measurementsOnPlane_.size() > 0))
      return true;
    errorOut << "KalmanFitterInfo::checkConsistency(): plane is nullptr" << std::endl;
    retVal = false;
  }

  TVector3 oTest = plane->getO(); // see if the plane object is there
  oTest *= 47.11;

  // if more than one measurement, check if they have the same dimensionality
  if (measurementsOnPlane_.size() > 1) {
    int dim = measurementsOnPlane_[0]->getState().GetNrows();
    for (unsigned int i = 1; i<measurementsOnPlane_.size(); ++i) {
      if(measurementsOnPlane_[i]->getState().GetNrows() != dim) {
        errorOut << "KalmanFitterInfo::checkConsistency(): measurementsOnPlane_ do not all have the same dimensionality" << std::endl;
        retVal = false;
      }
    }
    if (dim == 0) {
      errorOut << "KalmanFitterInfo::checkConsistency(): measurementsOnPlane_ have dimension 0" << std::endl;
      retVal = false;
    }
  }

  // see if everything else is defined wrt this plane and rep_
  int dim = rep_->getDim(); // check dim
  if (referenceState_) {
    if(referenceState_->getPlane() != plane) {
      errorOut << "KalmanFitterInfo::checkConsistency(): referenceState_ is not defined with the correct plane " << referenceState_->getPlane().get() << " vs. " << plane.get() << std::endl;
      retVal = false;
    }
    if (referenceState_->getRep() != rep_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): referenceState_ is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (referenceState_->getState().GetNrows() != dim) {
      errorOut << "KalmanFitterInfo::checkConsistency(): referenceState_ does not have the right dimension!" << std::endl;
      retVal = false;
    }
  }

  if (forwardPrediction_) {
    if(forwardPrediction_->getPlane() != plane) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ is not defined with the correct plane" << std::endl;
      retVal = false;
    }
    if(forwardPrediction_->getRep() != rep_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (forwardPrediction_->getState().GetNrows() != dim || forwardPrediction_->getCov().GetNrows() != dim) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ does not have the right dimension!" << std::endl;
      retVal = false;
    }
  }
  if (forwardUpdate_) {
    if(forwardUpdate_->getPlane() != plane) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ is not defined with the correct plane" << std::endl;
      retVal = false;
    }
    if(forwardUpdate_->getRep() != rep_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (forwardUpdate_->getState().GetNrows() != dim || forwardUpdate_->getCov().GetNrows() != dim) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ does not have the right dimension!" << std::endl;
      retVal = false;
    }
  }

  if (backwardPrediction_) {
    if(backwardPrediction_->getPlane() != plane) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ is not defined with the correct plane" << std::endl;
      retVal = false;
    }
    if(backwardPrediction_->getRep() != rep_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (backwardPrediction_->getState().GetNrows() != dim || backwardPrediction_->getCov().GetNrows() != dim) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ does not have the right dimension!" << std::endl;
      retVal = false;
    }
  }
  if (backwardUpdate_) {
    if(backwardUpdate_->getPlane() != plane) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ is not defined with the correct plane" << std::endl;
      retVal = false;
    }
    if(backwardUpdate_->getRep() != rep_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if (backwardUpdate_->getState().GetNrows() != dim || backwardUpdate_->getCov().GetNrows() != dim) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ does not have the right dimension!" << std::endl;
      retVal = false;
    }
  }

  for (std::vector<MeasurementOnPlane*>::const_iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    if((*it)->getPlane() != plane) {
      errorOut << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct plane" << std::endl;
      retVal = false;
    }
    if((*it)->getRep() != rep_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct TrackRep" << std::endl;
      retVal = false;
    }
    if ((*it)->getState().GetNrows() == 0) {
      errorOut << "KalmanFitterInfo::checkConsistency(): measurement has dimension 0!" << std::endl;
      retVal = false;
    }
  }

  if (flags == nullptr or !flags->hasFlags("U")) { // if predictions have not been pruned
    // see if there is an update w/o prediction or measurement
    if (forwardUpdate_ && !forwardPrediction_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o forwardPrediction_" << std::endl;
      retVal = false;
    }


    if (backwardUpdate_ && !backwardPrediction_) {
      errorOut << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o backwardPrediction_" << std::endl;
      retVal = false;
    }

    if (flags == nullptr or !flags->hasFlags("M")) {
      if (forwardUpdate_ && measurementsOnPlane_.size() == 0) {
        errorOut << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o measurement" << std::endl;
        retVal = false;
      }

      if (backwardUpdate_ && measurementsOnPlane_.size() == 0) {
        errorOut << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o measurement" << std::endl;
        retVal = false;
      }
    }
  }


  return retVal;
}


// Modified from auto-generated Streamer to correctly deal with smart pointers.
void KalmanFitterInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::KalmanFitterInfo.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::KalmanFitterInfo thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsFitterInfo baseClass0;
      baseClass0::Streamer(R__b);
      int flag;
      R__b >> flag;
      deleteForwardInfo();
      deleteBackwardInfo();
      deleteReferenceInfo();
      deleteMeasurementInfo();
      if (flag & 1) {
        referenceState_.reset(new ReferenceStateOnPlane());
        referenceState_->Streamer(R__b);
        referenceState_->setPlane(getPlane());
        // rep needs to be fixed up
      }
      if (flag & (1 << 1)) {
        forwardPrediction_.reset(new MeasuredStateOnPlane());
        forwardPrediction_->Streamer(R__b);
        forwardPrediction_->setPlane(getPlane());
        // rep needs to be fixed up
      }
      if (flag & (1 << 2)) {
        forwardUpdate_.reset(new KalmanFittedStateOnPlane());
        forwardUpdate_->Streamer(R__b);
        forwardUpdate_->setPlane(getPlane());
        // rep needs to be fixed up
      }
      if (flag & (1 << 3)) {	
        backwardPrediction_.reset(new MeasuredStateOnPlane());
        backwardPrediction_->Streamer(R__b);
        backwardPrediction_->setPlane(getPlane());
        // rep needs to be fixed up
      }
      if (flag & (1 << 4)) {	
        backwardUpdate_.reset(new KalmanFittedStateOnPlane());
        backwardUpdate_->Streamer(R__b);
        backwardUpdate_->setPlane(getPlane());
        // rep needs to be fixed up
      }
      {
        std::vector<genfit::MeasurementOnPlane*,std::allocator<genfit::MeasurementOnPlane*> > &R__stl =  measurementsOnPlane_;
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::MeasurementOnPlane));
        if (R__tcl1==0) {
          Error("measurementsOnPlane_ streamer","Missing the TClass object for genfit::MeasurementOnPlane!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        R__stl.reserve(R__n);
        for (R__i = 0; R__i < R__n; R__i++) {
          genfit::MeasurementOnPlane* R__t = new MeasurementOnPlane();
          R__t->Streamer(R__b);
          R__t->setPlane(getPlane());
          R__stl.push_back(R__t);
        }
      }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
     R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
     //This works around a msvc bug and should be harmless on other platforms
     typedef genfit::AbsFitterInfo baseClass0;
     baseClass0::Streamer(R__b);
     // "!!" forces the value to 1 or 0 (pointer != 0 or pointer == 0),
     // this value is then written as a bitfield.
     int flag = ((!!referenceState_)
		 | (!!forwardPrediction_ << 1)
		 | (!!forwardUpdate_ << 2)
		 | (!!backwardPrediction_ << 3)
		 | (!!backwardUpdate_ << 4));
     R__b << flag;
     if (flag & 1)
       referenceState_->Streamer(R__b);
     if (flag & (1 << 1))
       forwardPrediction_->Streamer(R__b);
     if (flag & (1 << 2))
       forwardUpdate_->Streamer(R__b);
     if (flag & (1 << 3))
       backwardPrediction_->Streamer(R__b);
     if (flag & (1 << 4))
       backwardUpdate_->Streamer(R__b);
     {
       std::vector<genfit::MeasurementOnPlane*,std::allocator<genfit::MeasurementOnPlane*> > &R__stl =  measurementsOnPlane_;
       int R__n=(&R__stl) ? int(R__stl.size()) : 0;
       R__b << R__n;
       if(R__n) {
         std::vector<genfit::MeasurementOnPlane*,std::allocator<genfit::MeasurementOnPlane*> >::iterator R__k;
         for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
           (*R__k)->Streamer(R__b);
         }
       }
     }
     R__b.SetByteCount(R__c, kTRUE);
  }
}


} /* End of namespace genfit */
