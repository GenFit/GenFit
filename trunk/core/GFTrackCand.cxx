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
#include "GFException.h"
#include "TDatabasePDG.h"

#include <algorithm>
#include <iostream>
#include <utility>

ClassImp(GFTrackCand)

GFTrackCand::GFTrackCand() : 
  fMcTrackId(-1),
  fPdg(0),
  fState6D(6),
  fCov6D(-1.0*TMatrixDSym(TMatrixDSym::kUnit,TMatrixDSym(6))),
  fQ(0) 
{
  ;
}

GFTrackCand::~GFTrackCand() {
  for (unsigned int i=0; i<fHits.size(); ++i) {
    delete fHits[i];
  }
  fHits.clear();
}


GFTrackCand::GFTrackCand( const GFTrackCand& other ) :
  fMcTrackId(other.fMcTrackId),
  fPdg(other.fPdg),
  fState6D(other.fState6D),
  fCov6D(other.fCov6D),
  fQ(other.fQ)
{
  // deep copy
  fHits.reserve(other.fHits.size());
  for (unsigned int i=0; i<other.fHits.size(); ++i) {
    fHits.push_back( new GFTrackCandHit(*(other.fHits[i])) );
  }
}

GFTrackCand&
GFTrackCand::operator=( const GFTrackCand& other ){
  fMcTrackId = other.fMcTrackId;
  fPdg = other.fPdg;
  fState6D = other.fState6D;
  fCov6D = other.fCov6D;
  fQ = other.fQ;

  for (unsigned int i=0; i<fHits.size(); ++i) {
    delete fHits[i];
  }
  fHits.clear();
  fHits.reserve(other.fHits.size());
  for (unsigned int i=0; i<other.fHits.size(); ++i) {
    fHits.push_back( new GFTrackCandHit(*(other.fHits[i])) );
  }

  return *this;
}


void 
GFTrackCand::addHit(int detId, int hitId, int planeId, double rho)
{
	fHits.push_back(new GFTrackCandHit(detId, hitId, planeId, rho));
}

std::vector<int>
GFTrackCand::getHitIDs(int detId) const {
  std::vector<int> result;
  for(unsigned int i=0; i<fHits.size(); ++i){
    if(detId==-2 || fHits[i]->getDetId() == detId) {
      result.push_back(fHits[i]->getHitId());
    }
  }
  return result;
}

std::vector<int>
GFTrackCand::getDetIDs() const {
  std::vector<int> result;
  for(unsigned int i=0; i<fHits.size(); ++i){
    result.push_back(fHits[i]->getDetId());
  }
  return result;
}

std::vector<double>
GFTrackCand::getRhos() const {
  std::vector<double> result;
  for(unsigned int i=0; i<fHits.size(); ++i){
    result.push_back(fHits[i]->getRho());
  }
  return result;
}

std::set<int>
GFTrackCand::getUniqueDetIDs() const {
  std::set<int> retVal;
  for (unsigned int i = 0; i < fHits.size(); ++i) {
    retVal.insert(fHits[i]->getDetId());
  }
  return retVal;
}


void
GFTrackCand::reset()
{
  for (unsigned int i=0; i<fHits.size(); ++i) {
    delete fHits[i];
  }
  fHits.clear();
}


bool GFTrackCand::hitInTrack(int detId, int hitId) const
{
	for (unsigned int i = 0; i < fHits.size(); ++i){
		if (detId == fHits[i]->getDetId() && hitId == fHits[i]->getHitId())
		  return true;
	}
	return false;	
}


bool operator== (const GFTrackCand& lhs, const GFTrackCand& rhs){
	if(lhs.getNHits() != rhs.getNHits()) return false;
	for (unsigned int i = 0; i < lhs.getNHits(); ++i){
	  if (lhs.getHit(i) != rhs.getHit(i)) return false;
	}
	return true;
}


void GFTrackCand::Print(const Option_t* option) const {
	std::cout << "======== GFTrackCand::print ========\n";
	std::cout << "mcTrackId=" << fMcTrackId << "\n";
	std::cout << "seed values for 6D state and cov: " << std::endl;
	fState6D.Print(option);
	fCov6D.Print(option);
	std::cout << "q" << fQ << "\n";
	std::cout << "PDG code= " << fPdg << "\n";
  for(unsigned int i=0; i<fHits.size(); ++i){
    fHits[i]->Print();
	}
}


void GFTrackCand::append(const GFTrackCand& rhs){
	for(unsigned int i=0; i<rhs.getNHits(); ++i){
		addHit(rhs.getHit(i));
	}
}


void GFTrackCand::sortHits(){
	std::stable_sort(fHits.begin(), fHits.end(), compareTrackCandHits);
}


void GFTrackCand::sortHits(const std::vector<unsigned int>& indices){

	const unsigned int nHits(getNHits());
	if (indices.size() != nHits){
		abort();
		GFException exc("GFTrackCand::sortHits ==> Size of indices != number of hits!",__LINE__,__FILE__);
		throw exc;
	}

	//these containers will hold the sorted results. They are created to avoid probably slower in-place sorting
	std::vector<GFTrackCandHit*> sortedHits(nHits);
	for (unsigned int i=0; i<nHits; ++i){
		sortedHits[i] = fHits[indices[i]];
	}
	//write the changes back to the private data members:
	fHits = sortedHits;
}
