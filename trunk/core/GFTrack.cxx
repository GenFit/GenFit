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
#include "GFException.h"
#include "GFAbsRecoHitComparator.h"
#include "TVirtualGeoTrack.h"

GFTrack::GFTrack(GFAbsTrackRep* defaultRep, bool smooth) 
  : fTrackReps(NULL),fCardinal_rep(0), fNextHitToFit(0), fSmooth(smooth)
{
  addTrackRep(defaultRep);
}

GFTrack::GFTrack() 
  : fTrackReps(NULL), fCardinal_rep(0), fNextHitToFit(0), fSmooth(false)
{
  ;
}

GFTrack::~GFTrack() {
  if(fTrackReps!=NULL){
    for(unsigned int i=0;i<getNumReps();i++) {
      delete fTrackReps->At(i);
    }
    delete fTrackReps;
  }
  for(unsigned int i=0;i<fHits.size();i++) {
    if(fHits[i]!=NULL) delete fHits[i];
  }
  for(unsigned int i=0;i<fBookkeeping.size();++i){
    if(fBookkeeping[i]!=NULL) delete fBookkeeping[i];
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
    assert(_tr.fBookkeeping[i]!= NULL) ;
    fBookkeeping.push_back(new GFBookkeeping(*(_tr.fBookkeeping[i])));
  }
  fRepAtHit = _tr.fRepAtHit;
}

GFTrack& GFTrack::operator=(const GFTrack& _tr) {
  if (this == &_tr)
    return *this;
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
    if(fBookkeeping[i]!=NULL) delete fBookkeeping[i];
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
    assert(_tr.fBookkeeping[i]!= NULL) ;
    fBookkeeping.push_back(new GFBookkeeping(*(_tr.fBookkeeping[i])));
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
    if(fBookkeeping[i]!=NULL) delete fBookkeeping[i];
  }
  for(unsigned int i=0;i<fHits.size();i++) {
    if(fHits[i]!=NULL) delete fHits[i];
  }
  fHits.clear();
  fRepAtHit.clear();
  fBookkeeping.clear();
}

void
GFTrack::mergeHits(GFTrack* trk){
  unsigned int nhits=trk->getNumHits();
  for(unsigned int i=0;i<nhits;++i){
    int detId, hitId;
    trk->getCand().getHit(i,detId,hitId);
    GFAbsRecoHit* hit=trk->getHit(i);
    addHit(hit,detId,hitId);
  }
  trk->fHits.clear();
}


void GFTrack::sortHits(){

  unsigned int nHits(getNumHits());

  std::vector< std::pair<unsigned int, GFAbsRecoHit*> > pv;
  pv.reserve(nHits);

  for (unsigned int i=0; i<nHits; ++i) {
      pv.push_back( std::pair<unsigned int, GFAbsRecoHit*>(i, fHits[i]) ) ;
  }

  std::stable_sort(pv.begin(), pv.end(), GFAbsRecoHitComparator());

  // get the indices -> now we know at which position which hit is
  std::vector<unsigned int> indices;
  indices.reserve(nHits);
  for (unsigned int i=0; i<nHits; ++i){
    indices.push_back(pv[i].first);
  }

  // now we have to do the actual sorting of everything

  // sort fHits
  std::vector<GFAbsRecoHit*> sortedHits;
  sortedHits.reserve(nHits);
  for (unsigned int i=0; i<nHits; ++i){
    sortedHits.push_back(fHits[indices[i]]);
  }
  fHits = sortedHits;

  // sort trackCand
  if (nHits == fCand.getNHits()){
    fCand.sortHits(indices);
  }
  else {
    GFException exc("GFTrack::sortHits ==> Cannot sort GFTrackCand accordingly since it has not the same number of hits as the GFTrack.",__LINE__,__FILE__);
    throw exc;
  }

  // sorting the bookkeeping is not yet supported (and probably wouldn't make sense either)
  clearBookkeeping();

  // update other member variables
  for (unsigned int i=0; i<fRepAtHit.size(); ++i){
    fRepAtHit[i] = std::find(indices.begin(), indices.end(), fRepAtHit[i]) - indices.begin();
  }

  // reset fNextHitToFit
  fNextHitToFit = 0;

}


void
GFTrack::setCandidate(const GFTrackCand& cand, bool doreset)
{
  fCand = cand;
  // reset fits
  if(doreset) {
    for(unsigned int i=0;i<getNumReps();i++) {
      (reinterpret_cast<GFAbsTrackRep*>(fTrackReps->At(i)))->reset();
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
    const GFDetPlane& pl( fHits[i]->getDetPlane(rep) );
    TVector3 pos = rep->getPos(pl);
    std::cout<<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<std::endl;
    geotrk->AddPoint(pos.X(),pos.Y(),pos.Z(),0);
  }// end loop over hits
}


void 
GFTrack::getResiduals(int detId, // which detector?
		    unsigned int dim,   // which projection?
		    unsigned int repid,   // which trackrep ?
		    std::vector<double>& result) const
{

  unsigned int nhits=getNumHits();
  if(repid>=getNumReps())return;
  GFAbsTrackRep* rep=getTrackRep(repid);//->clone();
  assert(rep->getState()==getTrackRep(repid)->getState());
  for(unsigned int ih=0; ih<nhits; ++ih){// loop over hits
    int anid, dummy;
    fCand.getHit(ih,anid,dummy); // check if this is a hit we want to look at
    if(anid==detId){
      GFAbsRecoHit* hit=getHit(ih);
      // extrapolate trackrep there
      int repDim=rep->getDim();
      TVectorD state(repDim);
      TMatrixDSym cov(repDim);
      const GFDetPlane& pl( hit->getDetPlane(rep) );
      
      rep->extrapolate(pl,state,cov);

      const TMatrixD& H( hit->getHMatrix(rep) );
      TVectorD m;
      TMatrixDSym V;
      hit->getMeasurement(rep,pl,state,cov,m,V);
      double res=(m-(H*state))[dim];

      result.push_back(res);
    } 
  }
}


void GFTrack::printBookkeeping(){
  std::cout << "GFTrack::printBookkeeping()" << std::endl;
  for(unsigned int i=0;i<getNumReps();++i){
    std::cout << "trackRep " << i << ":" << std::endl;    
    fBookkeeping[i]->Print();
  }

}

void GFTrack::Print(const Option_t* option) const{
  for(unsigned int i=0;i<getNumReps();++i){
    std::cout << "TrackRep " << i << " (defined at hit " << getRepAtHit(i) << "):\n";
    getTrackRep(i)->Print(option);
    fBookkeeping[i]->Print(option);
  }
  std::cout << "GFTrack has " << getNumHits() << " detector hits. ";
  if (fSmooth) std::cout << "Smoothing is enabled.";
  else std::cout << "Smoothing is disabled.";
  std::cout << std::endl;
}


void GFTrack::getHitMap(std::map<GFAbsRecoHit*, unsigned int>& hitMap) const {
  hitMap.clear();
  unsigned int nHits = getNumHits();

  for (unsigned int i=0; i<nHits; ++i){
    hitMap.insert(std::pair<GFAbsRecoHit*, unsigned int>(fHits[i], i));
  }

  return;
}


bool GFTrack::getHitsByPlane(std::vector< std::vector<int> >& retVal){

  retVal.clear();

  //this method can only be called when all hits have been loaded
  unsigned int nHits = fCand.getNHits();
  assert(fHits.size() == nHits);
  if(fHits.size() <= 1) return false;

  int detId, hitId, lastDetId, planeId, lastPlaneId;

  // get info for first hit
  fCand.getHitWithPlane(0,detId,hitId,planeId);

  lastPlaneId = planeId;
  lastDetId = detId;
  std::vector<int> vec;
  retVal.push_back(vec);
  retVal[0].reserve(10);
  retVal[0].push_back(0);

  // loop over hits
  for(unsigned int i=1; i<nHits; ++i){
    fCand.getHitWithPlane(i,detId,hitId,planeId);

    // if the next hit is in the same detector and has the same plane id (but not the default plane id -1), group them together
    if(planeId != -1 && planeId == lastPlaneId && detId == lastDetId){
      retVal.back().push_back(i);
    }
    else{
      lastPlaneId = planeId;
      lastDetId = detId;
      retVal.push_back(vec);
      retVal.back().reserve(10);
      retVal.back().push_back(i);
    }
  }

  return true;
}


void
GFTrack::blowUpCovs(double blowUpFactor){
  int nreps=getNumReps();
  for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=getTrackRep(irep);

    //dont do it for already compromsied reps, since they wont be fitted anyway
    if(arep->getStatusFlag()==0) {
      TMatrixDSym cov = arep->getCov();
      cov *= blowUpFactor;
      arep->setCov(cov);
    }
  }
}


ClassImp(GFTrack)



