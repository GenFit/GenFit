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
#include "GFFastKalman.h"

#include "assert.h"
#include <iostream>
#include <sstream>

#include "TMath.h"

#include "../Eigen/Dense"
#include "../Eigen/LU"

#include "GFTrack.h"
#include "GFAbsRecoHit.h"
#include "GFAbsTrackRep.h"
#include "GFException.h"
#include "GFTools.h"

#define COVEXC "cov_is_zero"
//#define DEBUG


void checkSym(TMatrixT<double> eval, int line){

	for(int i = 1; i < eval.GetNrows(); i++) for(int j = i; j < eval.GetNcols(); j++)
		if((eval(i,j) - eval(j,i)) >= 10e-10){
			GFException exc("COV is not symmetrical",line,__FILE__);
			throw exc;
		}
}



#define MAXROW 7																				// maximal dimensions of the Eigen Matrices. Choose wisely!
#define MAXCOL 7																				// Choose them big enough to fulfill the dimensional needs of all the matrices used within the Kalman Fitter, but small enough to be as fast as possible.
																								// Calculations are slowed down extremely when these are chosen as too big.

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 1, MAXROW, MAXCOL> EigMat;		// typedef for easier use of Eigen Matrices ('Eigen::Dynamic' lol...)

/*TMatrixT<double> calcGain(const TMatrixT<double>& cov,
						const TMatrixT<double>& HitCov,
						const TMatrixT<double>& H);
*/

void calcGain(const EigMat& Ecov, const EigMat& EHitCov, const EigMat& EH, EigMat& gain);

  /** @brief this returns the reduced chi2 increment for a hit
   */
/*
  double chi2Increment(const TMatrixT<double>& r,const TMatrixT<double>& H,
		       const TMatrixT<double>& cov,const TMatrixT<double>& V);
*/

  double chi2Increment(const EigMat& Er, const EigMat& EH, const EigMat& Ecov, const EigMat& EV);


  //void copyMatrices(EigMat& eMat, const TMatrixT<double>& TMat, const int& rows, const int& cols);


GFFastKalman::GFFastKalman():fInitialDirection(1),fNumIt(3),fBlowUpFactor(500.){



}


GFFastKalman::~GFFastKalman(){;}

void GFFastKalman::processTrack(GFTrack* trk){
#ifdef DEBUG
        std::cout<<"GFFastKalman::processTrack "<<std::endl;
#endif

  fSmooth = trk->getSmoothing();

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
      trk->getBK(i)->bookMatrices("fUpSt");
      trk->getBK(i)->bookMatrices("fUpCov");
      trk->getBK(i)->bookMatrices("bUpSt");
      trk->getBK(i)->bookMatrices("bUpCov");
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
GFFastKalman::switchDirection(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){
    trk->getTrackRep(i)->switchDirection();
  }
}

void GFFastKalman::blowUpCovs(GFTrack* trk){
  trk->blowUpCovs(fBlowUpFactor);
}

void
GFFastKalman::fittingPass(GFTrack* trk, int direction){
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
    //H_temp.ResizeTo(trk->getHit(ihit)->getHMatrix(arep).GetNrows(), trk->getHit(ihit)->getHMatrix(arep).GetNcols());
    if(arep->getStatusFlag()==0) {
      try {
#ifdef DEBUG
        std::cout<<"++++++++++++++++++++++++++++++++++++++++\n";
        std::cout<<"GFFastKalman::fittingPass - process rep nr. "<<irep<<" and hit nr. "<<ihit<<std::endl;
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
/*
double GFFastKalman::chi2Increment(const TMatrixT<double>& r,const TMatrixT<double>& H,
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
*/

double
GFFastKalman::getChi2Hit(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  // get prototypes for matrices
  int repDim=rep->getDim();

  EigMat state(repDim,1);											// This is the trick: double-declaration of matrices as Eigen and as ROOT objects.
  EigMat cov(repDim,repDim);										// But both work on the same data array!

						// Important: Use data array of Eigen matrices, because this is declared on the stack instead of on the heap
  						// Temporary dynamic matrices for use in TrackRep (needs to be resizeable)

  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	// which is important for Eigen to do fast calculations.
  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	// Array is unified.
  TMatrixT<double> state_temp;      // Temporary dynamic matrices for use in TrackRep (need to be resizeable)
  TMatrixT<double> cov_temp;

  GFDetPlane pl=hit->getDetPlane(rep);
  rep->extrapolate(pl,state_temp,cov_temp);


  memcpy(state.data(),state_temp.GetMatrixArray(),state.size() * sizeof(double));
  memcpy(cov.data(),cov_temp.GetMatrixArray(),cov.size() * sizeof(double));

  TMatrixT<double> H_temp = hit->getHMatrix(rep);
  TMatrixT<double> m_temp(state_temp.GetNcols(), H_temp.GetNcols());
  TMatrixT<double> V_temp(H_temp.GetNcols(), H_temp.GetNcols());

  EigMat H (H_temp.GetNrows(), H_temp.GetNcols());
  EigMat m (state_temp.GetNcols(), H_temp.GetNcols());
  EigMat V (V_temp.GetNcols(), V_temp.GetNcols());


  hit->getMeasurement(rep,pl,state_temp,cov_temp,m_temp,V_temp);

  memcpy(state.data(),state_temp.GetMatrixArray(),state.size() * sizeof(double));
  memcpy(cov.data(),cov_temp.GetMatrixArray(),cov.size() * sizeof(double));
  memcpy(m.data(),m_temp.GetMatrixArray(),m.size() * sizeof(double));
  memcpy(V.data(),V_temp.GetMatrixArray(),V.size() * sizeof(double));
  memcpy(H.data(),H_temp.GetMatrixArray(),H.size() * sizeof(double));

  EigMat res(m_temp.GetNrows(), m_temp.GetNcols());
  res = m-(H*state);
  //assert(Tm.GetNrows()>0);

  //this is where chi2 is calculated
  double chi2 = chi2Increment(res, H, cov, V);

  return chi2/m_temp.GetNrows();
}

  
void
GFFastKalman::processHit(GFTrack* tr, int ihit, int irep,int direction){
  GFAbsRecoHit* hit = tr->getHit(ihit);
  GFAbsTrackRep* rep = tr->getTrackRep(irep);

  int repDim=rep->getDim();

  EigMat state(repDim, 1);
  EigMat cov(repDim, repDim);

  state_temp.Use(repDim, 1, state.data());
  cov_temp.Use(repDim, repDim, cov.data());

  GFDetPlane pl;

  /* do an extrapolation, if the trackrep irep is not given
   * at this ihit position. This will usually be the case, but
   * not if the fit turnes around
   */
  if(ihit!=tr->getRepAtHit(irep)){
    // get the (virtual) detector plane
    pl=hit->getDetPlane(rep);
    //do the extrapolation
    rep->extrapolate(pl, state_ext, cov_ext);
    memcpy(state_temp.GetMatrixArray(),state_ext.GetMatrixArray(),state.size() * sizeof(double));
    memcpy(cov_temp.GetMatrixArray(),cov_ext.GetMatrixArray(),cov.size() * sizeof(double));
  }
  else{
    pl = rep->getReferencePlane();
    state_temp = rep->getState();
    cov_temp = rep->getCov();
  }


  if(cov(0,0)<1.E-50){ // diagonal elements must be >=0
    GFException exc(COVEXC,__LINE__,__FILE__);
    throw exc;
  }
  
#ifdef DEBUG
  std::cerr<<"GFFastKalman::processHit - state and cov prediction "<<std::endl;
  Tstate.Print();
  Tcov.Print();
#endif

  TMatrixT<double> H_temp = hit->getHMatrix(rep);
  TMatrixT<double> m_temp(state_temp.GetNcols(), H_temp.GetNcols());
  TMatrixT<double> V_temp(H_temp.GetNcols(), H_temp.GetNcols());

  EigMat H (H_temp.GetNrows(), H_temp.GetNcols());
  H_temp.Use(H_temp.GetNrows(), H_temp.GetNcols(), H.data());

  H_temp = hit->getHMatrix(rep);

  hit->getMeasurement(rep, pl, state_temp, cov_temp, m_temp, V_temp);

  EigMat m (H_temp.GetNrows(), 1);
  EigMat V (H_temp.GetNrows(), H_temp.GetNrows());

  memcpy(V.data(),V_temp.GetMatrixArray(),V.size() * sizeof(double));
 // m_temp.ResizeTo(H_temp.GetNrows(), 1);
  memcpy(m.data(),m_temp.GetMatrixArray(),m.size() * sizeof(double));

  EigMat res(m_temp.GetNrows(), m_temp.GetNcols());

  res = m-(H*state);

  // calculate kalman gain ------------------------------
  EigMat Gain ( cov_temp.GetNrows() , V_temp.GetNcols());

  calcGain(cov, V, H, Gain);

  //std::cout<<"Gain:\n" << Gain << std::endl;


  // calculate update -----------------------------------
  EigMat update(state_temp.GetNrows(), state_temp.GetNcols());


  update=Gain*res;


#ifdef DEBUG
  std::cout<<"Gain\n" << Gain << std::endl;
  std::cout<<"residual vector\n" << res << std::endl;
  std::cout<<"update = Gain*res\n" << update << std::endl;
#endif

  state+=update; // prediction overwritten!
  cov-=Gain*(H*cov);

  if(fSmooth) {
    if(direction == 1) {
    tr->getBK(irep)->setMatrix("fUpSt",ihit,state_temp);
    tr->getBK(irep)->setMatrix("fUpCov",ihit,cov_temp);
    if(rep->hasAuxInfo()) tr->getBK(irep)->setMatrix("fAuxInfo",ihit,*(rep->getAuxInfo(pl)));
    tr->getBK(irep)->setDetPlane("fPl",ihit,pl);
  } else {
    tr->getBK(irep)->setMatrix("bUpSt",ihit,state_temp);
    tr->getBK(irep)->setMatrix("bUpCov",ihit,cov_temp);
    if(rep->hasAuxInfo()) tr->getBK(irep)->setMatrix("bAuxInfo",ihit,*(rep->getAuxInfo(pl)));
    tr->getBK(irep)->setDetPlane("bPl",ihit,pl);
  }
  }


  res = m-(H*state);

  double chi2 = chi2Increment(res, H, cov, V);


  int ndf = m_temp.GetNrows();
  if (direction == -1) {
    rep->addChiSqu( chi2 );
  }
  if (direction == 1) {
    rep->addForwardChiSqu( chi2 );
  }
  rep->addNDF( ndf );


  rep->setData(state_temp,pl,&cov_temp);
  tr->setRepAtHit(irep,ihit);

#ifdef DEBUG
   std::cout<<"GFFastKalman::processHit - updated state and cov "<<std::endl;
   rep->getState().Print();
   rep->getCov().Print();
#endif


}


void calcGain(const EigMat& Ecov, const EigMat& EHitCov, const EigMat& EH, EigMat& gain){

	  gain = Ecov*(EH.transpose()*((EHitCov + EH*Ecov*(EH.transpose())).inverse()));

}


double chi2Increment(const EigMat& Er, const EigMat& EH, const EigMat& Ecov, const EigMat& EV){

	EigMat chisq (1,1);
	chisq = (Er.transpose())*((EV - EH*Ecov*(EH.transpose())).inverse())*Er;

	assert(chisq.size() == 1);
	if(TMath::IsNaN(chisq(0,0))){
	GFException exc("chi2 is nan",__LINE__,__FILE__);
	exc.setFatal();

	throw exc;
	}

	return chisq(0,0);
}

