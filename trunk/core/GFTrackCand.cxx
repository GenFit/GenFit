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
#include "GFTrackCand.h"
#include "TDatabasePDG.h"
#include <algorithm>
#include <iostream>

ClassImp(GFTrackCand)

GFTrackCand::GFTrackCand():fCurv(0),fDip(0),fInv(false), fQoverpSeed(0.),fMcTrackId(-1),fPdg(0){}

GFTrackCand::~GFTrackCand(){}

GFTrackCand::GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs)
  : fCurv(curv), fDip(dip), fInv(inv),fQoverpSeed(0.), fMcTrackId(-1),fPdg(0)
{
  assert(detIDs.size()==hitIDs.size());
  int n = detIDs.size();
  for ( int i; i != n; i++){
	  TrackCandHit aNewHit = {detIDs[i], hitIDs[i], 0.0, 0};
	  fTrackCandHits.push_back(aNewHit);
  }
}
GFTrackCand::GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs,std::vector<double> rhos)
  : fCurv(curv), fDip(dip), fInv(inv),fQoverpSeed(0.), fMcTrackId(-1),fPdg(0)
{
  assert(detIDs.size()==hitIDs.size());
  assert(detIDs.size()==rhos.size());
  int n = detIDs.size();
  for ( int i; i != n; i++){
	  TrackCandHit aNewHit = {detIDs[i], hitIDs[i], 0.0, 0};
	  fTrackCandHits.push_back(aNewHit);
  }
}

void 
GFTrackCand::addHit(unsigned int detId, unsigned int hitId, double rho, unsigned int planeId)
{
  TrackCandHit aNewHit = {detId, hitId, rho, planeId};
  fTrackCandHits.push_back(aNewHit);
}

std::vector<unsigned int>
GFTrackCand::GetHitIDs(int detId) const {
	std::vector<unsigned int> result;
	int n=fTrackCandHits.size();
  if(detId<0){ // return hits from all detectors
	    for(int i = 0; i != n; ++i){
	    	  result.push_back(fTrackCandHits[i].fHitId);
	    }
  }
  else {
    for(int i = 0; i != n; ++i){
      if(fTrackCandHits[i].fDetId == unsigned(detId)){
    	  result.push_back(fTrackCandHits[i].fHitId);
      }
    }
  }
  return result;
}

std::vector<double>
GFTrackCand::GetRhos(int detId) const {
	std::vector<double> result;
	int n=fTrackCandHits.size();
  if(detId<0){ // return hits from all detectors
	    for(int i = 0; i != n; ++i){
	    	  result.push_back(fTrackCandHits[i].fRho);
	    }
  }
  else {
    for(int i = 0; i != n; ++i){
      if(fTrackCandHits[i].fDetId == unsigned(detId)){
    	  result.push_back(fTrackCandHits[i].fRho);
      }
    }
  }
  return result;
}

void
GFTrackCand::reset()
{
	fTrackCandHits.clear();
}

bool GFTrackCand::HitInTrack(unsigned int detId, unsigned int hitId) const
{
	int nHits = fTrackCandHits.size();
	for (int i = 0; i != nHits; i++){
		if (detId == fTrackCandHits[i].fDetId && hitId == fTrackCandHits[i].fHitId){
			return true;
		}
	}
	return false;	
}

bool operator== (const GFTrackCand& lhs, const GFTrackCand& rhs){
  if(lhs.getNHits() != rhs.getNHits()){
	  return false;
  }
  int n = lhs.getNHits();
  for (int i = 0; i != n; ++i){
	  if ( lhs.fTrackCandHits[i].fDetId != rhs.fTrackCandHits[i].fDetId || lhs.fTrackCandHits[i].fHitId !=  rhs.fTrackCandHits[i].fHitId ){
		  return false;
	  }
  }
  return true;
}

void GFTrackCand::Print(const Option_t* option) const {
  std::cout << "======== GFTrackCand::print ========" << std::endl;
  if(fMcTrackId>=0) std::cout << "mcTrackId=" << fMcTrackId << std::endl;
  std::cout << "seed values for pos,direction, and q/p: " << std::endl;
  fPosSeed.Print(option);
  fDirSeed.Print(option);
  std::cout << "q/p=" << fQoverpSeed << std::endl;
  std::cout << "detId|hitId|rho ";
  int n = getNHits();
  for( int i = 0; i != n; ++i){
    std::cout << fTrackCandHits[i].fDetId << "|" << fTrackCandHits[i].fHitId
	      << "|" << fTrackCandHits[i].fRho << " ";
  }
  std::cout << std::endl;
}

void GFTrackCand::append(const GFTrackCand& rhs){
  unsigned int detId;
  unsigned int hitId;
  double rho;
  int n = rhs.getNHits();
  for ( int i=0; i != n ;++i){
    rhs.getHit(i,detId,hitId,rho);
    addHit(detId,hitId,rho);
  }
}

void GFTrackCand::setComplTrackSeed(const TVector3& pos,const TVector3& mom, const int pdgCode, TVector3 posError, TVector3 dirError){
  fPosSeed=pos;fDirSeed=mom; fPdg = pdgCode; fPosError = posError; fDirError = dirError;
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(fPdg);
  double charge = part->Charge()/(3.);
  fQoverpSeed=charge/mom.Mag();
}

void GFTrackCand::sortHits(){
	std::sort(fTrackCandHits.begin(), fTrackCandHits.end(), compareRho );
}

bool GFTrackCand::compareRho( const TrackCandHit& lhs, const TrackCandHit& rhs){
  	return lhs.fRho < rhs.fRho;
}
