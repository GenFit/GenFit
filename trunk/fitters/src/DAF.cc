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
  setProbCut(0.01);
}

DAF::DAF(AbsKalmanFitter* kalman, double deltaPval, double deltaWeight)
  : AbsKalmanFitter(10, deltaPval), deltaWeight_(deltaWeight)
{
  kalman_.reset(kalman);
  kalman_->setMultipleMeasurementHandling(weightedAverage); // DAF makes no sense otherwise
  kalman_->setMaxIterations(1);

  if (dynamic_cast<KalmanFitterRefTrack*>(kalman_.get()) != NULL) {
    static_cast<KalmanFitterRefTrack*>(kalman_.get())->setRefitAll();
  }

  setAnnealingScheme(100, 0.1, 5); // also sets maxIterations_
  setProbCut(0.01);
}


void DAF::processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits) {

  if (debugLvl_ > 0) {
    std::cout<<"DAF::processTrack //////////////////////////////////////////////////////////////// \n";
  }

  KalmanFitStatus* status = 0;
  bool oneLastIter = false;

  double lastPval = -1;

  for(unsigned int iBeta = 0;; ++iBeta) {

    if (debugLvl_ > 0) {
      std::cout<<"DAF::processTrack, trackRep  " << rep << ", iteration " << iBeta+1 << ", beta = " << betas_.at(iBeta) << "\n";
    }

    kalman_->processTrackWithRep(tr, rep, resortHits);

    status = static_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
    status->setIsFittedWithDaf();
    status->setNumIterations(iBeta+1);


    // check break conditions

    if (! status->isFitted()){
      if (debugLvl_ > 0) {
        std::cout << "DAF::Kalman could not fit!\n";
      }
      status->setIsFitted(false);
      break;
    }

    if( oneLastIter == true){
      if (debugLvl_ > 0) {
        std::cout << "DAF::break after one last iteration\n";
      }
      status->setIsFitConvergedFully(status->getNFailedPoints() == 0);
      status->setIsFitConvergedPartially();
      break;
    }

    if(iBeta >= maxIterations_-1){
      status->setIsFitConvergedFully(false);
      status->setIsFitConvergedPartially(false);
      if (debugLvl_ > 0) {
        std::cout << "DAF::number of max iterations reached!\n";
      }
      break;
    }


    // get and update weights
    bool converged(false);
    try{
      converged = calcWeights(tr, rep, betas_.at(iBeta));
      if (!converged && iBeta >= minIterations_-1 &&
          status->getBackwardPVal() != 0 && abs(lastPval - status->getBackwardPVal()) < this->deltaPval_) {
        if (debugLvl_ > 0) {
          std::cout << "converged by Pval = " << status->getBackwardPVal() << " even though weights changed at iBeta = " << iBeta << std::endl;
        }
        converged = true;
      }
      lastPval = status->getBackwardPVal();
    } catch(Exception& e) {
      std::cerr<<e.what();
      e.info();
      //std::cerr << "calc weights failed" << std::endl;
      //mini_trk->getTrackRep(0)->setStatusFlag(1);
      status->setIsFitted(false);
      status->setIsFitConvergedFully(false);
      status->setIsFitConvergedPartially(false);
      break;
    }

    // check if converged
    if (iBeta >= minIterations_-1 && converged) {
      if (debugLvl_ > 0) {
        std::cout << "DAF::convergence reached in iteration " << iBeta+1 << " -> Do one last iteration with updated weights.\n";
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
  for ( int i = 1; i != 6; ++i){
    addProbCut(prob_cut, i);
  }
}

void DAF::addProbCut(const double prob_cut, const int measDim){
  if ( prob_cut > 1.0 || prob_cut < 0.0){
    Exception exc("DAF::addProbCut prob_cut is not between 0 and 1",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  if ( measDim < 1){
    Exception exc("DAF::addProbCut measDim must be > 0",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  chi2Cuts_[measDim] = ROOT::Math::chisquared_quantile_c( prob_cut, measDim);
}


void DAF::setBetas(double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8,double b9,double b10){
  betas_.clear();
  assert(b1>0);betas_.push_back(b1);
  if(b2>0){
    assert(b2<=b1);betas_.push_back(b2);
    if(b3>=0.) {
      assert(b3<=b2);betas_.push_back(b3);
      if(b4>=0.) {
        assert(b4<=b3);betas_.push_back(b4);
        if(b5>=0.) {
          assert(b5<=b4);betas_.push_back(b5);
          if(b6>=0.) {
            assert(b6<=b5);betas_.push_back(b6);
            if(b7>=0.) {
              assert(b7<=b6);betas_.push_back(b7);
              if(b8>=0.) {
                assert(b8<=b7);betas_.push_back(b8);
                if(b9>=0.) {
                  assert(b9<=b8);betas_.push_back(b9);
                  if(b10>=0.) {
                    assert(b10<=b9);betas_.push_back(b10);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  minIterations_ = betas_.size();
  maxIterations_ = betas_.size() + 4;
  betas_.resize(maxIterations_,betas_.back()); //make sure main loop has a maximum of maxIterations_ and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.
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
    std::cout<< betas_.at(i) << ", ";
  }*/
}


bool DAF::calcWeights(Track* tr, const AbsTrackRep* rep, double beta) {

  if (debugLvl_ > 0) {
    std::cout<<"DAF::calcWeights \n";
  }

  bool converged(true);
  double maxAbsChange(0);

  const std::vector< TrackPoint* >& trackPoints = tr->getPointsWithMeasurement();
  for (std::vector< TrackPoint* >::const_iterator tp = trackPoints.begin(); tp != trackPoints.end(); ++tp) {
    if (! (*tp)->hasFitterInfo(rep)) {
      continue;
    }
    AbsFitterInfo* fi = (*tp)->getFitterInfo(rep);
    if (dynamic_cast<KalmanFitterInfo*>(fi) == NULL){
      Exception exc("DAF::getWeights ==> can only use KalmanFitterInfos",__LINE__,__FILE__);
      throw exc;
    }
    KalmanFitterInfo* kfi = static_cast<KalmanFitterInfo*>(fi);

    if (kfi->areWeightsFixed()) {
      if (debugLvl_ > 0) {
        std::cout<<"weights are fixed, continue \n";
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
          std::cout<<"chi2 = " << chi2 << "\n";
        }

        double norm = 1./sqrt(twoPiN * detV);

        phi[j] = norm*exp(-0.5*chi2/beta);
        phi_sum += phi[j];
        //std::cerr << "hitDim " << hitDim << " fchi2Cuts[hitDim] " << fchi2Cuts[hitDim] << std::endl;
        double cutVal = chi2Cuts_[hitDim];
        assert(cutVal>1.E-6);
        //the following assumes that in the competing hits could have different V otherwise calculation could be simplified
        phi_cut += norm*exp(-0.5*cutVal/beta);
      }
      catch(Exception& e) {
        std::cerr << e.what();
        e.info();
      }
    }

    for(unsigned int j=0; j<nMeas; j++) {
      double weight = phi[j]/(phi_sum+phi_cut);

      // check convergence
      double absChange(fabs(weight - kfi->getMeasurementOnPlane(j)->getWeight()));
      if (converged && absChange > deltaWeight_) {
        converged = false;
        if (absChange > maxAbsChange)
          maxAbsChange = absChange;
      }

      if (debugLvl_ > 0) {
        if (debugLvl_ > 1 || absChange > deltaWeight_) {
          std::cout<<"\t old weight: " << kfi->getMeasurementOnPlane(j)->getWeight();
          std::cout<<"\t new weight: " << weight;
        }
      }

      kfi->getMeasurementOnPlane(j)->setWeight(weight);
    }
  }

  if (debugLvl_ > 0) {
    std::cout << "\t  ";
    std::cout << "max abs weight change = " << maxAbsChange << "\n";
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
      {
         std::map<int,double> &R__stl =  chi2Cuts_;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            double R__t2;
            R__b >> R__t2;
            typedef int Value_t;
            std::pair<Value_t const, double > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
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
         int R__n=(&R__stl) ? int(R__stl.size()) : 0;
         R__b << R__n;
         if(R__n) {
            std::vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         std::map<int,double> &R__stl =  chi2Cuts_;
         int R__n=(&R__stl) ? int(R__stl.size()) : 0;
         R__b << R__n;
         if(R__n) {
            std::map<int,double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << ((*R__k).first );
            R__b << ((*R__k).second);
            }
         }
      }
      R__b << kalman_.get();
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
