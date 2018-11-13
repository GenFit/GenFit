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

#include "DAF.h"
#include "Exception.h"
#include "IO.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include "KalmanFitStatus.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

#include <assert.h>
#include <cmath>

//root stuff
#include <TBuffer.h>
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>



namespace genfit {

DAF::DAF(bool useRefKalman, double deltaPval, double deltaWeight)
  : AbsKalmanFitter(10, deltaPval), deltaWeight_(deltaWeight)
{
  if (useRefKalman) {
    kalman_.reset(new KalmanFitterRefTrack());
    static_cast<KalmanFitterRefTrack*>(kalman_.get())->setRefitAll();
  }
  else
    kalman_.reset(new KalmanFitter());

  kalman_->setMultipleMeasurementHandling(weightedAverage);
  kalman_->setMaxIterations(1);

  setAnnealingScheme(100, 0.1, 5); // also sets maxIterations_
  setProbCut(0.001);
}

DAF::DAF(AbsKalmanFitter* kalman, double deltaPval, double deltaWeight)
  : AbsKalmanFitter(10, deltaPval), deltaWeight_(deltaWeight)
{
  kalman_.reset(kalman);
  kalman_->setMultipleMeasurementHandling(weightedAverage); // DAF makes no sense otherwise
  kalman_->setMaxIterations(1);

  if (dynamic_cast<KalmanFitterRefTrack*>(kalman_.get()) != nullptr) {
    static_cast<KalmanFitterRefTrack*>(kalman_.get())->setRefitAll();
  }

  setAnnealingScheme(100, 0.1, 5); // also sets maxIterations_
  setProbCut(0.01);
}


void DAF::processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits) {

  if (debugLvl_ > 0) {
    debugOut<<"DAF::processTrack //////////////////////////////////////////////////////////////// \n";
  }

  KalmanFitStatus* status = 0;
  bool oneLastIter = false;

  double lastPval = -1;

  for(unsigned int iBeta = 0;; ++iBeta) {

    if (debugLvl_ > 0) {
      debugOut<<"DAF::processTrack, trackRep  " << rep << ", iteration " << iBeta+1 << ", beta = " << betas_.at(iBeta) << "\n";
    }

    kalman_->processTrackWithRep(tr, rep, resortHits);

    status = static_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
    status->setIsFittedWithDaf();
    status->setNumIterations(iBeta+1);


    // check break conditions

    if (! status->isFitted()){
      if (debugLvl_ > 0) {
        debugOut << "DAF::Kalman could not fit!\n";
      }
      status->setIsFitted(false);
      break;
    }

    if( oneLastIter == true){
      if (debugLvl_ > 0) {
        debugOut << "DAF::break after one last iteration\n";
      }
      status->setIsFitConvergedFully(status->getNFailedPoints() == 0);
      status->setIsFitConvergedPartially();
      break;
    }

    if(iBeta >= maxIterations_-1){
      status->setIsFitConvergedFully(false);
      status->setIsFitConvergedPartially(false);
      if (debugLvl_ > 0) {
        debugOut << "DAF::number of max iterations reached!\n";
      }
      break;
    }


    // get and update weights
    bool converged(false);
    try{
      converged = calcWeights(tr, rep, betas_.at(iBeta));
      if (!converged && iBeta >= minIterations_-1 &&
          status->getBackwardPVal() != 0 && fabs(lastPval - status->getBackwardPVal()) < this->deltaPval_) {
        if (debugLvl_ > 0) {
          debugOut << "converged by Pval = " << status->getBackwardPVal() << " even though weights changed at iBeta = " << iBeta << std::endl;
        }
        converged = true;
      }
      lastPval = status->getBackwardPVal();
    } catch(Exception& e) {
      errorOut<<e.what();
      e.info();
      //errorOut << "calc weights failed" << std::endl;
      //mini_trk->getTrackRep(0)->setStatusFlag(1);
      status->setIsFitted(false);
      status->setIsFitConvergedFully(false);
      status->setIsFitConvergedPartially(false);
      break;
    }

    // check if converged
    if (iBeta >= minIterations_-1 && converged) {
      if (debugLvl_ > 0) {
        debugOut << "DAF::convergence reached in iteration " << iBeta+1 << " -> Do one last iteration with updated weights.\n";
      }
      oneLastIter = true;
      status->setIsFitConvergedFully(status->getNFailedPoints() == 0);
      status->setIsFitConvergedPartially();
    }

  } // end loop over betas


  if (status->getForwardPVal() == 0. &&
      status->getBackwardPVal() == 0.) {
    status->setIsFitConvergedFully(false);
    status->setIsFitConvergedPartially(false);
  }

}


void DAF::setProbCut(const double prob_cut){
  for ( int i = 1; i < 7; ++i){
    addProbCut(prob_cut, i);
  }
}

void DAF::addProbCut(const double prob_cut, const int measDim){
  if ( prob_cut > 1.0 || prob_cut < 0.0){
    Exception exc("DAF::addProbCut prob_cut is not between 0 and 1",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  if ( measDim < 1 || measDim > 6 ){
    Exception exc("DAF::addProbCut measDim must be in the range [1,6]",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  chi2Cuts_[measDim] = ROOT::Math::chisquared_quantile_c( prob_cut, measDim);
}

void DAF::setAnnealingScheme(double bStart, double bFinal, unsigned int nSteps) {
  // The betas are calculated as a geometric series that takes nSteps
  // steps to go from bStart to bFinal.
  assert(bStart > bFinal);
  assert(bFinal > 1.E-10);
  assert(nSteps > 1);

  minIterations_ = nSteps;
  maxIterations_ = nSteps + 4;

  betas_.clear();

  for (unsigned int i=0; i<nSteps; ++i) {
    betas_.push_back(bStart * pow(bFinal / bStart, i / (nSteps - 1.)));
  }

  betas_.resize(maxIterations_,betas_.back()); //make sure main loop has a maximum of 10 iterations and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.

  /*for (unsigned int i=0; i<betas_.size(); ++i) {
    debugOut<< betas_.at(i) << ", ";
  }*/
}


bool DAF::calcWeights(Track* tr, const AbsTrackRep* rep, double beta) {

  if (debugLvl_ > 0) {
    debugOut<<"DAF::calcWeights \n";
  }

  bool converged(true);
  double maxAbsChange(0);

  const std::vector< TrackPoint* >& trackPoints = tr->getPointsWithMeasurement();
  for (std::vector< TrackPoint* >::const_iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {
    if (! (*tp)->hasFitterInfo(rep)) {
      continue;
    }
    AbsFitterInfo* fi = (*tp)->getFitterInfo(rep);
    if (dynamic_cast<KalmanFitterInfo*>(fi) == nullptr){
      Exception exc("DAF::getWeights ==> can only use KalmanFitterInfos",__LINE__,__FILE__);
      throw exc;
    }
    KalmanFitterInfo* kfi = static_cast<KalmanFitterInfo*>(fi);

    if (kfi->areWeightsFixed()) {
      if (debugLvl_ > 0) {
        debugOut<<"weights are fixed, continue \n";
      }
      continue;
    }

    unsigned int nMeas = kfi->getNumMeasurements();

    std::vector<double> phi(nMeas, 0.);
    double phi_sum = 0;
    double phi_cut = 0;
    for(unsigned int j=0; j<nMeas; j++) {

      try{
        const MeasurementOnPlane& residual = kfi->getResidual(j, true, true);
        const TVectorD& resid(residual.getState());
        TMatrixDSym Vinv(residual.getCov());
        double detV;
        tools::invertMatrix(Vinv, &detV); // can throw an Exception
        int hitDim = resid.GetNrows();
        // Needed for normalization, special cases for the two common cases,
        // shouldn't matter, but the original code made some efforts to make
        // this calculation faster, and it's not complex ...
        double twoPiN = 2.*M_PI;
        if (hitDim == 2)
          twoPiN *= twoPiN;
        else if (hitDim > 2)
          twoPiN = pow(twoPiN, hitDim);

        double chi2 = Vinv.Similarity(resid);
        if (debugLvl_ > 1) {
          debugOut<<"chi2 = " << chi2 << "\n";
        }

	// The common factor beta is eliminated.
        double norm = 1./sqrt(twoPiN * detV);

        phi[j] = norm*exp(-0.5*chi2/beta);
	phi_sum += phi[j];
        //errorOut << "hitDim " << hitDim << " fchi2Cuts[hitDim] " << fchi2Cuts[hitDim] << std::endl;
        double cutVal = chi2Cuts_[hitDim];
        assert(cutVal>1.E-6);
        //the following assumes that in the competing hits could have different V otherwise calculation could be simplified
        phi_cut += norm*exp(-0.5*cutVal/beta);
      }
      catch(Exception& e) {
        errorOut << e.what();
        e.info();
      }
    }

    for(unsigned int j=0; j<nMeas; j++) {
      double weight = phi[j]/(phi_sum+phi_cut);
      //debugOut << phi_sum << " " << phi_cut << " " << weight << std::endl;

      // check convergence
      double absChange(fabs(weight - kfi->getMeasurementOnPlane(j)->getWeight()));
      if (converged && absChange > deltaWeight_) {
        converged = false;
        if (absChange > maxAbsChange)
          maxAbsChange = absChange;
      }

      if (debugLvl_ > 0) {
        if (debugLvl_ > 1 || absChange > deltaWeight_) {
          debugOut<<"\t old weight: " << kfi->getMeasurementOnPlane(j)->getWeight();
          debugOut<<"\t new weight: " << weight;
        }
      }

      kfi->getMeasurementOnPlane(j)->setWeight(weight);
    }
  }

  if (debugLvl_ > 0) {
    debugOut << "\t  ";
    debugOut << "max abs weight change = " << maxAbsChange << "\n";
  }

  return converged;
}


// Customized from generated Streamer.
void DAF::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::DAF.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::DAF thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsKalmanFitter baseClass0;
      baseClass0::Streamer(R__b);
      R__b >> deltaWeight_;
      // weights_ are only of intermediate use -> not saved
      {
         std::vector<double> &R__stl =  betas_;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      if (R__v == 1) {
 	 // Old versions kept the chi2Cuts_ in a map.
         // We ignore non-sensical dimensionalities when reading it again.
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
	    memset(chi2Cuts_, 0, sizeof(chi2Cuts_));
            int R__t;
            R__b >> R__t;
            double R__t2;
            R__b >> R__t2;
	    if (R__t >= 1 && R__t <= 6)
	      chi2Cuts_[R__t] = R__t2;
	 }
      } else {
	char n_chi2Cuts;  // should be six
	R__b >> n_chi2Cuts;
	assert(n_chi2Cuts == 6);  // Cannot be different as long as sanity prevails.
	chi2Cuts_[0] = 0; // nonsensical.
	R__b.ReadFastArray(&chi2Cuts_[1], n_chi2Cuts);
      }
      AbsKalmanFitter *p;
      R__b >> p;
      kalman_.reset(p);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsKalmanFitter baseClass0;
      baseClass0::Streamer(R__b);
      R__b << deltaWeight_;
      {
         std::vector<double> &R__stl =  betas_;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            std::vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
	R__b << (char)6;  // Number of chi2Cuts_
	R__b.WriteFastArray(&chi2Cuts_[1], 6);
      }
      R__b << kalman_.get();
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
