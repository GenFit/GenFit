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

#include "assert.h"
#include <iostream>
#include <sstream>

#include "TMath.h"

#include "GFTrack.h"
#include "GFAbsRecoHit.h"
#include "GFAbsTrackRep.h"
#include "GFException.h"
#include "GFTools.h"

#define COVEXC "cov_is_zero"
//#define DEBUG

GFKalman::GFKalman():fInitialDirection(1),fNumIt(3),fBlowUpFactor(500.){;}

GFKalman::~GFKalman(){;}

void GFKalman::processTrack(GFTrack* trk){
#ifdef DEBUG
        std::cout<<"GFKalman::processTrack "<<std::endl;
#endif

  fSmooth = trk->getSmoothing();
  fSmoothFast = trk->getSmoothingFast();

  int nreps = trk->getNumReps();
  for(int i=0; i<nreps; i++) {
    trk->getBK(i)->setNhits(trk->getNumHits());
    if(fSmooth) {
      std::vector<std::string> mat_keys = trk->getBK(i)->getMatrixKeys();
      bool already_there = false;
      for(unsigned int j=0; j<mat_keys.size(); j++) {
        if(mat_keys.at(j) == "fUpSt") already_there = true;
      }
      if(already_there) continue;
      trk->getBK(i)->bookNumbers("fExtLen"); // extrapolated length from last hit in forward direction
      trk->getBK(i)->bookMatrices("fUpSt");
      trk->getBK(i)->bookMatrices("fUpCov");
      trk->getBK(i)->bookNumbers("bExtLen"); // extrapolated length from last hit in backward direction
      trk->getBK(i)->bookMatrices("bUpSt");
      trk->getBK(i)->bookMatrices("bUpCov");
      if(fSmoothFast) {
          trk->getBK(i)->bookMatrices("fSt");
          trk->getBK(i)->bookMatrices("fCov");
          trk->getBK(i)->bookMatrices("bSt");
          trk->getBK(i)->bookMatrices("bCov");
      }
      trk->getBK(i)->bookGFDetPlanes("fPl");
      trk->getBK(i)->bookGFDetPlanes("bPl");
      if(trk->getTrackRep(i)->hasAuxInfo()) {
        trk->getBK(i)->bookMatrices("fAuxInfo");
        trk->getBK(i)->bookMatrices("bAuxInfo");
      }
    }
  }

  int direction=fInitialDirection;
  assert(direction==1 || direction==-1);
  //  trk->clearGFBookkeeping();
  trk->clearRepAtHit();

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
    
    //save first and last plane,state&cov after the fitting pass
    if(direction==1){//forward at last hit
      for(int i=0; i<nreps; ++i){
        trk->getTrackRep(i)->
          setLastPlane( trk->getTrackRep(i)->getReferencePlane() );
        trk->getTrackRep(i)->
          setLastState( trk->getTrackRep(i)->getState() );
        trk->getTrackRep(i)->
          setLastCov( trk->getTrackRep(i)->getCov() );
      }
    }
    else{//backward at first hit
      for(int i=0; i<nreps; ++i){
        trk->getTrackRep(i)->
          setFirstPlane( trk->getTrackRep(i)->getReferencePlane() );
        trk->getTrackRep(i)->
          setFirstState( trk->getTrackRep(i)->getState() );
        trk->getTrackRep(i)->
          setFirstCov( trk->getTrackRep(i)->getCov() );
      }
    }

    //switch direction of fitting and also inside all the reps
    if(direction==1) direction=-1;
    else direction=1;
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

void GFKalman::blowUpCovs(GFTrack* trk){
  trk->blowUpCovs(fBlowUpFactor);
}

void
GFKalman::fittingPass(GFTrack* trk, int direction){
  //loop over hits
  unsigned int nhits=trk->getNumHits();
  unsigned int starthit=trk->getNextHitToFit();

  int nreps=trk->getNumReps();
  int ihit=(int)starthit;
  
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
    //    GFAbsRecoHit* ahit=trk->getHit(ihit);
    // loop over reps
    for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    if(arep->getStatusFlag()==0) {
      try {
#ifdef DEBUG
        std::cout<<"++++++++++++++++++++++++++++++++++++++++\n";
        std::cout<<"GFKalman::fittingPass - process rep nr. "<<irep<<" and hit nr. "<<ihit<<std::endl;
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
  trk->setNextHitToFit(ihit-2*direction);
  //trk->printGFBookkeeping();
}

double GFKalman::chi2Increment(const TMatrixT<double>& r,const TMatrixT<double>& H,
           const TMatrixT<double>& cov,const TMatrixT<double>& V){

  // residuals covariances:R=(V - HCH^T)
  TMatrixT<double> R(V);
  TMatrixT<double> covsum1(cov,TMatrixT<double>::kMultTranspose,H);
  TMatrixT<double> covsum(H,TMatrixT<double>::kMult,covsum1);

  R-=covsum;

  // chisq= r^TR^(-1)r
  TMatrixT<double> Rinv;
  GFTools::invertMatrix(R,Rinv);
  
  
  TMatrixT<double> residTranspose(r);
  residTranspose.T();
  TMatrixT<double> chisq=residTranspose*(Rinv*r);
  assert(chisq.GetNoElements()==1);

  if(TMath::IsNaN(chisq[0][0])){
    GFException exc("chi2 is nan",__LINE__,__FILE__);
    exc.setFatal();
    std::vector< TMatrixT<double> > matrices;
    matrices.push_back(r);
    matrices.push_back(V);
    matrices.push_back(R);
    matrices.push_back(cov);
    exc.setMatrices("r, V, R, cov",matrices);
    throw exc;
  }

  return chisq[0][0];
}


double
GFKalman::getChi2Hit(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  // get prototypes for matrices
  int repDim=rep->getDim();
  TMatrixT<double> state(repDim,1);
  TMatrixT<double> cov(repDim,repDim);;
  GFDetPlane pl=hit->getDetPlane(rep);
  rep->extrapolate(pl,state,cov);


  TMatrixT<double> H = hit->getHMatrix(rep);
  TMatrixT<double> m,V;
  hit->getMeasurement(rep,pl,state,cov,m,V);

  TMatrixT<double> res = m-(H*state);
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
  TMatrixT<double> state(repDim,1);
  TMatrixT<double> cov(repDim,repDim);;
  GFDetPlane pl;

  double extLen(0.);

  /* do an extrapolation, if the trackrep irep is not given
   * at this ihit position. This will usually be the case, but
   * not if the fit turnes around
   */
  if(ihit!=tr->getRepAtHit(irep)){
    // get the (virtual) detector plane
    pl=hit->getDetPlane(rep);
    //do the extrapolation
    extLen = rep->extrapolate(pl,state,cov);
  }
  else{
    pl = rep->getReferencePlane();
    state = rep->getState();
    cov = rep->getCov();
    extLen = 0.;
  }
  
  if(cov[0][0]<1.E-50){ // diagonal elements must be >=0
    GFException exc(COVEXC,__LINE__,__FILE__);
    throw exc;
  }

  if(fSmooth && fSmoothFast) {
    if(direction == 1) {
	    tr->getBK(irep)->setMatrix("fSt",ihit,state);
	    tr->getBK(irep)->setMatrix("fCov",ihit,cov);
	    if(rep->hasAuxInfo()) tr->getBK(irep)->setMatrix("fAuxInfo",ihit,*(rep->getAuxInfo(pl)));
	    tr->getBK(irep)->setDetPlane("fPl",ihit,pl);
	  } else {
	    tr->getBK(irep)->setMatrix("bSt",ihit,state);
	    tr->getBK(irep)->setMatrix("bCov",ihit,cov);
	    if(rep->hasAuxInfo()) tr->getBK(irep)->setMatrix("bAuxInfo",ihit,*(rep->getAuxInfo(pl)));
	    tr->getBK(irep)->setDetPlane("bPl",ihit,pl);
	  }
  }
  
#ifdef DEBUG
  std::cerr<<"GFKalman::processHit - state and cov prediction "<<std::endl;
  state.Print();
  cov.Print();
#endif

  TMatrixT<double> H(hit->getHMatrix(rep));
  TMatrixT<double> m, V;
  hit->getMeasurement(rep,pl,state,cov,m,V);
  TMatrixT<double> res = m-(H*state);

  // calculate kalman gain ------------------------------
  TMatrixT<double> Gain(calcGain(cov,V,H));

  // calculate update -----------------------------------
  TMatrixT<double> update=Gain*res;

#ifdef DEBUG
  std::cout<<"Gain"; Gain.Print();
  std::cout<<"residual vector"; res.Print();
  std::cout<<"update = Gain*res"; update.Print();
#endif

  state+=update; // prediction overwritten!
  cov-=Gain*(H*cov);

  if(fSmooth) {
    if(direction == 1) {
      tr->getBK(irep)->setNumber("fExtLen",ihit,extLen);
      tr->getBK(irep)->setMatrix("fUpSt",ihit,state);
      tr->getBK(irep)->setMatrix("fUpCov",ihit,cov);
      if(rep->hasAuxInfo()) tr->getBK(irep)->setMatrix("fAuxInfo",ihit,*(rep->getAuxInfo(pl)));
      tr->getBK(irep)->setDetPlane("fPl",ihit,pl);
    } else {
	    tr->getBK(irep)->setNumber("bExtLen",ihit,extLen);
      tr->getBK(irep)->setMatrix("bUpSt",ihit,state);
      tr->getBK(irep)->setMatrix("bUpCov",ihit,cov);
      if(rep->hasAuxInfo()) tr->getBK(irep)->setMatrix("bAuxInfo",ihit,*(rep->getAuxInfo(pl)));
      tr->getBK(irep)->setDetPlane("bPl",ihit,pl);
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

  /*
  if(direction==1){
    tr->getBK(irep)->setNumber("fChi2",ihit,chi2/ndf);
  }
  else{
    tr->getBK(irep)->setNumber("bChi2",ihit,chi2/ndf);
  }
  */

  // if we survive until here: update TrackRep
  //rep->setState(state);
  //rep->setCov(cov);
  //rep->setReferencePlane(pl);

  rep->setData(state,pl,&cov);
  tr->setRepAtHit(irep,ihit);

#ifdef DEBUG
   std::cout<<"GFKalman::processHit - updated state and cov "<<std::endl;
   rep->getState().Print();
   rep->getCov().Print();
#endif
}


TMatrixT<double>
GFKalman::calcGain(const TMatrixT<double>& cov, 
           const TMatrixT<double>& HitCov,
           const TMatrixT<double>& H){

  // calculate covsum (V + HCH^T)
  TMatrixT<double> covsum1(cov,TMatrixT<double>::kMultTranspose,H);
  TMatrixT<double> covsum(H,TMatrixT<double>::kMult,covsum1);

  covsum+=HitCov;
  
  // invert
  TMatrixT<double> covsumInv;
  GFTools::invertMatrix(covsum,covsumInv);
  
  // calculate gain
  TMatrixT<double> gain1(H,TMatrixT<double>::kTransposeMult,covsumInv);
  TMatrixT<double> gain(cov,TMatrixT<double>::kMult,gain1);

  return gain;
}







