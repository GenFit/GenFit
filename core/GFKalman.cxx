/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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
#include "GFKalman.h"

#include <assert.h>
#include <sstream>

#include <TMath.h>

#include "GFTrack.h"
#include "GFBookkeeping.h"
#include "GFException.h"
#include "GFTools.h"

#define COVEXC "cov_is_zero"
//#define DEBUG

GFKalman::GFKalman()
  : fInitialDirection(1), fNumIt(3)
{
  ;
}

GFKalman::~GFKalman(){;}

void GFKalman::processTrack(GFTrack* trk){
#ifdef DEBUG
  std::cout<<"GFKalman::processTrack with " << fNumIt << " iterations." <<std::endl;
#endif

  initBookkeeping(trk);
  trk->clearRepAtHit();

  int direction=fInitialDirection;
  assert(direction==1 || direction==-1);

  int nreps = trk->getNumReps();

  /*why is there a factor of two here (in the for statement)?
    Because we consider one full iteration to be one back and
    one forth fitting pass */
  for(int ipass=0; ipass<2*fNumIt; ipass++){
    if(ipass>0) blowUpCovs(trk);

    // reset X/X0 before last fitting pass
    if(ipass==(2*fNumIt)-1) {
      for(int i=0; i<nreps; ++i) {
        trk->getTrackRep(i)->resetXX0();
      }
    }


    if(direction==1){
      trk->setNextHitToFit(0);
    }
    else {
      trk->setNextHitToFit(trk->getNumHits()-1);
    }

    fittingPass(trk,direction);
    
    //save first and last plane, state and cov after the fitting pass
    if(direction==1){
      for(int i=0; i<nreps; ++i){
        GFAbsTrackRep* rep = trk->getTrackRep(i);
        rep->setLastPlane( rep->getReferencePlane() );
        rep->setLastState( rep->getState() );
        rep->setLastCov( rep->getCov() );
      }
    }
    else{
      for(int i=0; i<nreps; ++i){
        GFAbsTrackRep* rep = trk->getTrackRep(i);
        rep->setFirstPlane( rep->getReferencePlane() );
        rep->setFirstState( rep->getState() );
        rep->setFirstCov( rep->getCov() );
      }
    }

    //switch direction of fitting and also inside all the reps
    direction *= -1;
    switchDirection(trk);
  }
  
  return;
}

void
GFKalman::switchDirection(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){
    trk->getTrackRep(i)->switchDirection();
  }
}

void
GFKalman::fittingPass(GFTrack* trk, int direction){
  //loop over hits
  unsigned int nhits = trk->getNumHits();
  unsigned int starthit = trk->getNextHitToFit();

  int nreps = trk->getNumReps();
  int ihit = (int)starthit;


  for(int irep=0; irep<nreps; ++irep) {
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    if(arep->getStatusFlag()==0) {
      //clear chi2 sum and ndf sum in track reps
        if (direction == -1){
          arep->setChiSqu(0.);
        }
        if (direction == 1){
          arep->setForwardChiSqu(0.);
        }
      arep->setNDF(0);
      //clear failedHits and outliers
      trk->getBK(irep)->clearFailedHits();
    }
  }
  while((ihit<(int)nhits && direction==1) || (ihit>-1 && direction==-1)){
    // loop over reps
    for(int irep=0; irep<nreps; ++irep){
      GFAbsTrackRep* arep=trk->getTrackRep(irep);
      if(arep->getStatusFlag()==0) {
        try {
#ifdef DEBUG
          std::cout<<"++++++++++++++++++++++++++++++++++++++++\n";
          std::cout<<"GFKalman::fittingPass, direction = " << direction << ";  process rep nr. "<<irep<<" and hit nr. "<<ihit<<std::endl;
#endif
          processHit(trk,ihit,irep,direction);
        }
        catch(GFException& e) {
          trk->addFailedHit(irep,ihit);
          std::cerr << e.what();
          e.info();
          if(e.isFatal()) {
            arep->setStatusFlag(1);
            continue; // go to next rep immediately
          }
        }
      }
    }// end loop over reps

    ihit+=direction;
  }// end loop over hits
}

double GFKalman::chi2Increment(const TVectorD& r,
                               const TMatrixD& H,
                               const TMatrixDSym& cov,
                               const TMatrixDSym& V){

  // residuals covariances: Rinv = R^-1 = (V - HCH^T)^-1
  TMatrixDSym Rinv(cov);
  Rinv.Similarity(H);

  Rinv -= V; // this is now  -R = HCH^T - V
  GFTools::invertMatrix(Rinv); // this is now  -R^(-1)
  double chisq = Rinv.Similarity(r); // -chisq = r^T -R^(-1) r
  chisq *= -1.;

  if(TMath::IsNaN(chisq)){
    GFException exc("chi2 is nan",__LINE__,__FILE__);
    exc.setFatal();
    std::vector<TMatrixD> matrices;
    matrices.push_back(V);
    matrices.push_back(Rinv);
    matrices.push_back(cov);
    exc.setMatrices("V, R^(-1), cov",matrices);
    throw exc;
  }

  return chisq;
}


double
GFKalman::getChi2Hit(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  // get prototypes for matrices
  int repDim=rep->getDim();
  TVectorD state(repDim);
  TMatrixDSym cov(repDim);;
  const GFDetPlane& pl(hit->getDetPlane(rep));
  rep->extrapolate(pl,state,cov);


  const TMatrixD& H = hit->getHMatrix(rep);
  TVectorD m;
  TMatrixDSym V;
  hit->getMeasurement(rep,pl,state,cov,m,V);

  TVectorD res = m-(H*state);
  assert(res.GetNrows()>0);

  //this is where chi2 is calculated
  double chi2 = chi2Increment(res,H,cov,V);

  return chi2/res.GetNrows();
}

  
void
GFKalman::processHit(GFTrack* tr, int ihit, int irep,int direction){
  GFAbsRecoHit* hit = tr->getHit(ihit);
  GFAbsTrackRep* rep = tr->getTrackRep(irep);

  // get prototypes for matrices
  int repDim = rep->getDim();
  TVectorD state(repDim);
  TMatrixDSym cov(repDim);
  const GFDetPlane* ppl;

  double extLen(0.);

  /* do an extrapolation, if the trackrep irep is not given
   * at this ihit position. This will usually be the case, but
   * not if the fit turnes around
   */
  if(ihit!=tr->getRepAtHit(irep)){
    // get the (virtual) detector plane
    ppl=&hit->getDetPlane(rep);
    //do the extrapolation
    extLen = rep->extrapolate(*ppl,state,cov);
  }
  else{
    ppl = &rep->getReferencePlane();
    state = rep->getState();
    cov = rep->getCov();
    extLen = 0.;
  }
  const GFDetPlane& pl = *ppl;

  if(cov(0,0)<1.E-50){ // diagonal elements must be >=0
    GFException exc(COVEXC,__LINE__,__FILE__);
    throw exc;
  }

  GFBookkeeping* bk = tr->getBK(irep);
  if(tr->getSmoothing()) { // save predictions
    if(direction == 1) {
	    bk->setVector(GFBKKey_fSt,ihit,state);
	    bk->setSymMatrix(GFBKKey_fCov,ihit,cov);
	  } else {
	    bk->setVector(GFBKKey_bSt,ihit,state);
	    bk->setSymMatrix(GFBKKey_bCov,ihit,cov);
	  }
  }
  
#ifdef DEBUG
  std::cerr<<"GFKalman::processHit - state and cov prediction "<<std::endl;
  state.Print();
  cov.Print();
#endif

  const TMatrixD& H(hit->getHMatrix(rep));
  TVectorD m;
  TMatrixDSym V;
  hit->getMeasurement(rep,pl,state,cov,m,V);
  TVectorD res = m-(H*state);

  // calculate kalman gain ------------------------------
  // calculate covsum (V + HCH^T)
  TMatrixDSym HcovHt(cov);
  HcovHt.Similarity(H);
  
  // invert
  //TMatrixDSym covSum(V + HcovHt);
  // instead of:
  //  TMatrixDSym covSum(V, TMatrixDSym::kPlus, HcovHt);
  // because kPlus constructor doesn't work in root up to at least 5.34.
  // Bug report: <https://savannah.cern.ch/bugs/index.php?98605>
  TMatrixDSym covSumInv(V + HcovHt);
  GFTools::invertMatrix(covSumInv);

  TMatrixD CHt(cov, TMatrixD::kMultTranspose, H);
  TVectorD update = TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res;
#ifdef DEBUG
  std::cout<<"residual vector"; res.Print();
  std::cout<<"update = Gain*res"; update.Print();
#endif

  state+=update; // prediction overwritten!

  // And the new covariance matrix:
  covSumInv.Similarity(CHt);
  cov-=covSumInv;  // Cnew = C - C Ht (V + H C Ht)^-1 H C

  if(tr->getSmoothing()) {
    if(direction == 1) {
      bk->setNumber(GFBKKey_fExtLen,ihit,extLen);
      bk->setVector(GFBKKey_fUpSt,ihit,state);
      bk->setSymMatrix(GFBKKey_fUpCov,ihit,cov);
      bk->setDetPlane(GFBKKey_fPl,ihit,pl);
      if(rep->hasAuxInfo()) bk->setMatrix(GFBKKey_fAuxInfo,ihit,*(rep->getAuxInfo(pl)));
    } else {
	    bk->setNumber(GFBKKey_bExtLen,ihit,extLen);
      bk->setVector(GFBKKey_bUpSt,ihit,state);
      bk->setSymMatrix(GFBKKey_bUpCov,ihit,cov);
      bk->setDetPlane(GFBKKey_bPl,ihit,pl);
      if(rep->hasAuxInfo()) bk->setMatrix(GFBKKey_bAuxInfo,ihit,*(rep->getAuxInfo(pl)));
    }
  }

  // calculate filtered chisq
  // filtered residual
  res = m-(H*state);
  double chi2 = chi2Increment(res,H,cov,V);
  int ndf = res.GetNrows();
  if (direction == -1) {
    rep->addChiSqu( chi2 );
  }
  if (direction == 1) {
    rep->addForwardChiSqu( chi2 );
  }
  rep->addNDF( ndf );

  // update TrackRep
  rep->setData(state,pl,&cov);
  tr->setRepAtHit(irep,ihit);

#ifdef DEBUG
   std::cout<<"GFKalman::processHit - updated state and cov "<<std::endl;
   rep->getState().Print();
   rep->getCov().Print();
#endif
}


void
GFKalman::initBookkeeping(GFTrack* trk) const {

  int nreps = trk->getNumReps();
  for(int i=0; i<nreps; ++i) {
    GFBookkeeping* bk = trk->getBK(i);
    bk->setNhits(trk->getNumHits());
    if(trk->getSmoothing()) {
      const std::vector<GFBKKey>& keys = bk->getGFDetPlaneKeys();
      bool already_there = false;
      for(unsigned int j=0; j<keys.size(); ++j) {
        if(keys[j] == GFBKKey_bPl) {
          already_there = true;
          break;
        }
      }
      if(!already_there) {
        bk->bookVectors(GFBKKey_fSt);       // state prediction in forward direction
        bk->bookSymMatrices(GFBKKey_fCov);  // covariance prediction in forward direction
        bk->bookVectors(GFBKKey_fUpSt);     // updated state in forward direction
        bk->bookSymMatrices(GFBKKey_fUpCov);// updated covariance in forward direction
        bk->bookGFDetPlanes(GFBKKey_fPl);   // detector plane in forward direction
        bk->bookNumbers(GFBKKey_fExtLen);   // extrapolated length from last hit in forward direction

        bk->bookVectors(GFBKKey_bSt);       // state prediction in backward direction
        bk->bookSymMatrices(GFBKKey_bCov);  // covariance prediction in backward direction
        bk->bookVectors(GFBKKey_bUpSt);     // updated state in backward direction
        bk->bookSymMatrices(GFBKKey_bUpCov);// updated covariance in backward direction
        bk->bookGFDetPlanes(GFBKKey_bPl);   // detector plane in backward direction
        bk->bookNumbers(GFBKKey_bExtLen);   // extrapolated length from last hit in backward direction

        if(trk->getTrackRep(i)->hasAuxInfo()) {
          bk->bookMatrices(GFBKKey_fAuxInfo); // aux info in forward direction
          bk->bookMatrices(GFBKKey_bAuxInfo); // aux info in backward direction
        }
      }
    }// end if(trk->getSmoothing())
  } // end loop over reps
}







