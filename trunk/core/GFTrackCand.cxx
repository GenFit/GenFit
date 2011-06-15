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

#include <algorithm>
#include <iostream>

ClassImp(GFTrackCand)

GFTrackCand::GFTrackCand():fCurv(0),fDip(0),fInv(false), fQoverpSeed(0.),fMcTrackId(-1){}

GFTrackCand::~GFTrackCand(){}

GFTrackCand::GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs)
  : fDetId(detIDs),fHitId(hitIDs),fCurv(curv), fDip(dip), fInv(inv),fQoverpSeed(0.), fMcTrackId(-1)
{
  assert(fDetId.size()==fHitId.size());
  fRho.resize(fDetId.size(),0.);
}
GFTrackCand::GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs,std::vector<double> rhos)
  : fDetId(detIDs),fHitId(hitIDs),fRho(rhos),fCurv(curv), fDip(dip), fInv(inv),fQoverpSeed(0.), fMcTrackId(-1)
{
  assert(fDetId.size()==fHitId.size());
  assert(fDetId.size()==fRho.size());
}

void 
GFTrackCand::addHit(unsigned int detId, unsigned int hitId, double rho, unsigned int planeId)
{
  fDetId.push_back(detId);
  fHitId.push_back(hitId);
  fPlaneId.push_back(planeId);
  fRho.push_back(rho);
}

std::vector<unsigned int> 
GFTrackCand::GetHitIDs(int detId){
  if(detId<0){ // return hits from all detectors
    return fHitId;
  }
  else {
    std::vector<unsigned int> result;
    unsigned int n=fHitId.size();
    for(unsigned int i=0;i<n;++i){
      if(fDetId[i]==(unsigned int)detId)result.push_back(fHitId[i]);
    }
    return result;
  }
}

void
GFTrackCand::reset()
{
  fDetId.clear();fHitId.clear();
}

bool GFTrackCand::HitInTrack(unsigned int detId, unsigned int hitId)
{
	for (unsigned int i = 0; i < fDetId.size(); i++){
		if (detId == fDetId[i])
			if (hitId == fHitId[i])
				return true;
	}
	return false;	
}

bool operator== (const GFTrackCand& lhs, const GFTrackCand& rhs){
  if(lhs.getNHits()!=rhs.getNHits()) return false;
  bool result=std::equal(lhs.fDetId.begin(),lhs.fDetId.end(),rhs.fDetId.begin());
  result &=std::equal(lhs.fHitId.begin(),lhs.fHitId.end(),rhs.fHitId.begin());
  return result;
}

void GFTrackCand::Print() const {
  std::cout << "======== GFTrackCand::print ========" << std::endl;
  if(fMcTrackId>=0) std::cout << "mcTrackId=" << fMcTrackId << std::endl;
  std::cout << "seed values for pos,direction, and q/p: " << std::endl;
  fPosSeed.Print();
  fDirSeed.Print();
  std::cout << "q/p=" << fQoverpSeed << std::endl;
  assert(fDetId.size()==fHitId.size());
  std::cout << "detId|hitId|rho ";
  for(unsigned int i=0;i<fDetId.size();++i){
    std::cout << fDetId.at(i) << "|" << fHitId.at(i) 
	      << "|" << fRho.at(i) << " ";
  }
  std::cout << std::endl;
}

void GFTrackCand::append(const GFTrackCand& rhs){
  unsigned int detId,hitId;
  double rho;
  for(unsigned int i=0;i<rhs.getNHits();++i){
    rhs.getHit(i,detId,hitId,rho);
    addHit(detId,hitId,rho);
  }


}
