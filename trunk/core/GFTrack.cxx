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
#include <assert.h>
#include <iostream>

#include "GFTrack.h"
#include "TVirtualGeoTrack.h"

GFTrack::GFTrack(GFAbsTrackRep* defaultRep, bool smooth) 
  : fTrackReps(NULL),fCardinal_rep(0), fNextHitToFit(0), fSmooth(false)
{
  addTrackRep(defaultRep);
  fSmooth = smooth;
}

GFTrack::GFTrack() 
  : fTrackReps(NULL), fCardinal_rep(0), fNextHitToFit(0), fSmooth(false)
{
  //trackReps = new TObjArray(defNumTrackReps);
}

GFTrack::~GFTrack() {
  if(fTrackReps!=NULL){
    for(unsigned int i=0;i<getNumReps();i++) {
      delete fTrackReps->At(i);
    }
    delete fTrackReps;
  }
  for(unsigned int i=0;i<fHits.size();i++) {
    if(fHits.at(i)!=NULL) delete fHits.at(i);
  }
  for(unsigned int i=0;i<fBookkeeping.size();++i){
    if(fBookkeeping.at(i)!=NULL) delete fBookkeeping.at(i);
  }
}

GFTrack::GFTrack(const GFTrack& _tr) {
  fCand=_tr.fCand;
  fCardinal_rep=_tr.fCardinal_rep;
  fNextHitToFit=_tr.fNextHitToFit;
  fSmooth=_tr.fSmooth;
  for(unsigned int i=0;i<_tr.getNumHits();i++) {
    fHits.push_back((_tr.getHit(i))->clone());
  }
  fTrackReps = NULL;
  for(unsigned int i=0; i<_tr.getNumReps();i++) {
    addTrackRep( (_tr.getTrackRep(i))->clone() );
  }
  for(unsigned int i=0; i<fBookkeeping.size(); ++i) delete fBookkeeping[i];
  fBookkeeping.clear();

  for(unsigned int i=0;i<_tr.fBookkeeping.size();++i){
    assert(_tr.fBookkeeping.at(i)!= NULL) ;
    fBookkeeping.push_back(new GFBookkeeping(*(_tr.fBookkeeping.at(i))));
  }
  fRepAtHit = _tr.fRepAtHit;
}

GFTrack& GFTrack::operator=(const GFTrack& _tr) {
  if(fTrackReps!=NULL){
    for(unsigned int i=0;i<getNumReps();i++) {
      delete fTrackReps->At(i);
    }
    delete fTrackReps;
    fTrackReps=NULL;
  }
  for(unsigned int i=0;i<fHits.size();i++) {
    delete fHits[i];
  }
  for(unsigned int i=0;i<fBookkeeping.size();++i){
    if(fBookkeeping.at(i)!=NULL) delete fBookkeeping.at(i);
  }

  for(unsigned int i=0;i<_tr.getNumReps();++i){
    addTrackRep(_tr.getTrackRep(i)->clone());
  }
  fCand=_tr.fCand;
  fCardinal_rep=_tr.fCardinal_rep;
  fNextHitToFit=_tr.fNextHitToFit;
  fSmooth=_tr.fSmooth;
  for(unsigned int i=0;i<_tr.getNumHits();i++) {
    fHits.push_back((_tr.getHit(i))->clone());
  }

  //clear the empty bookeeping objs made by addTrackRep and copy the others
  for(unsigned int i=0; i<fBookkeeping.size(); ++i) delete fBookkeeping[i];
  fBookkeeping.clear();
  for(unsigned int i=0;i<_tr.fBookkeeping.size();++i){
    assert(_tr.fBookkeeping.at(i)!= NULL) ;
    fBookkeeping.push_back(new GFBookkeeping(*(_tr.fBookkeeping.at(i))));
  }
  fRepAtHit = _tr.fRepAtHit;


  return *this;
}


void
GFTrack::reset(){
  if(fTrackReps!=NULL){
    for(unsigned int i=0;i<getNumReps();i++) {
      if(fTrackReps->At(i)!=NULL) delete fTrackReps->At(i);
    }
  }
  for(unsigned int i=0;i<fBookkeeping.size();++i){
    if(fBookkeeping.at(i)!=NULL) delete fBookkeeping.at(i);
  }
  for(unsigned int i=0;i<fHits.size();i++) {
    if(fHits.at(i)!=NULL) delete fHits.at(i);
  }
  fHits.clear();
  fRepAtHit.clear();
  fBookkeeping.clear();
}

void
GFTrack::mergeHits(GFTrack* trk){
  unsigned int nhits=trk->getNumHits();
  for(unsigned int i=0;i<nhits;++i){
    unsigned int detId;
    unsigned int hitId;
    trk->getCand().getHit(i,detId,hitId);
    GFAbsRecoHit* hit=trk->getHit(i);
    addHit(hit,detId,hitId);
  }
  trk->fHits.clear();
}


void 
GFTrack::setCandidate(const GFTrackCand& cand, bool doreset)
{
  fCand=cand;
  // reset fits
  if(doreset) {
    for(unsigned int i=0;i<getNumReps();i++) {
      ((GFAbsTrackRep*)fTrackReps->At(i))->reset();
    }
  }
}

void 
GFTrack::fillGeoTrack(TVirtualGeoTrack* geotrk,unsigned int repid) const
{
  GFAbsTrackRep* rep=getTrackRep(repid);
  unsigned int n=fCand.getNHits();
  rep->getState().Print();
  for(unsigned int i=0; i<n; ++i){// loop over hits
    GFDetPlane pl=fHits[i]->getDetPlane(rep);
    TVector3 pos=rep->getPos(pl);
    std::cout<<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<std::endl;
    geotrk->AddPoint(pos.X(),pos.Y(),pos.Z(),0);
  }// end loop over hits
}


void 
GFTrack::getResiduals(unsigned int detId, // which detector?
		    unsigned int dim,   // which projection?
		    unsigned int repid,   // which trackrep ?
		    std::vector<double>& result)
{
  unsigned int nhits=getNumHits();
  if(repid>=getNumReps())return;
  GFAbsTrackRep* rep=getTrackRep(repid);//->clone();
  assert(rep->getState()==getTrackRep(repid)->getState());
  for(unsigned int ih=0; ih<nhits; ++ih){// loop over hits
    unsigned int anid;
    unsigned int dummy;
    fCand.getHit(ih,anid,dummy); // check if this is a hit we want to look at
    if(anid==detId){
      GFAbsRecoHit* hit=getHit(ih);
      // extrapolate trackrep there
      int repDim=rep->getDim();
      TMatrixT<double> state(repDim,1);
      GFDetPlane pl=hit->getDetPlane(rep);
      
      rep->extrapolate(pl,state);
      //rep->setState(state);
      //rep->setReferencePlane(pl);
      double res=hit->residualVector(rep,state,pl)[dim][0];

      //std::cout<<res<<std::endl;

      result.push_back(res);
    } 
  }
}


void GFTrack::printBookkeeping(){
  std::cout << "GFTrack::printBookkeeping()" << std::endl;
  for(unsigned int i=0;i<getNumReps();++i){
    std::cout << "trackRep " << i << ":" << std::endl;    
    fBookkeeping.at(i)->Print();
  }

}

void GFTrack::Print(const Option_t* option) const{
  for(unsigned int i=0;i<getNumReps();++i){
    getTrackRep(i)->Print(option);
    fBookkeeping.at(i)->Print(option);
  }
  std::cout << "GFTrack has " << getNumHits() << " detector hits." << std::endl;
  
}


void GFTrack::getHitsByPlane(std::vector<std::vector<int>*>& retVal){
  for(int i=0;retVal.size();++i){
    delete retVal.at(i);
  }
  retVal.clear();
  //this method can only be called when all hits have been loaded
  assert(fHits.size()==fCand.getNHits());
  assert(fHits.size()>1);
  unsigned int detId,hitId,planeId;
  fCand.getHitWithPlane(0,detId,hitId,planeId);
  //  std::cout << "$$$ " << 0 << " " << detId << " " << hitId << " " << planeId << std::endl;
  unsigned int lastPlane=planeId;
  retVal.push_back(new std::vector<int>);
  retVal.at(0)->push_back(0);
  for(unsigned int i=1;i<fCand.getNHits();++i){
    fCand.getHitWithPlane(i,detId,hitId,planeId);
    //std::cout << "$$$ " << i << " " << detId << " " << hitId << " " << planeId << std::endl;
    if(lastPlane==planeId){
      retVal.at(retVal.size()-1)->push_back(i);
    }
    else{
      lastPlane=planeId;
      retVal.push_back(new std::vector<int>);
      retVal.at(retVal.size()-1)->push_back(i);
    }
  }
}


void
GFTrack::blowUpCovs(double blowUpFactor){
  int nreps=getNumReps();
  for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=getTrackRep(irep);

    //dont do it for already compromsied reps, since they wont be fitted anyway
    if(arep->getStatusFlag()==0) {
      TMatrixT<double> cov = arep->getCov();
      for(int i=0;i<cov.GetNrows();++i){
        for(int j=0;j<cov.GetNcols();++j){
          //if(i!=j){//off diagonal
          //  cov[i][j]=0.;
          //}
          //else{//diagonal
          cov[i][j] = cov[i][j] * blowUpFactor;
          //}
        }
      }
      arep->setCov(cov);
    }
  }
}


ClassImp(GFTrack)



