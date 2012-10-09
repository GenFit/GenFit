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

GFTrackCand::GFTrackCand():fCurv(0),fDip(0), fQ(0),fMcTrackId(-1),fPdg(0), fState6D(TMatrixD(TMatrixD::kZero,TMatrixD(6,1))),fCov6D(-1.0*TMatrixD(TMatrixD::kUnit,TMatrixD(6,6))){}

GFTrackCand::~GFTrackCand(){} //no members are pointers

void 
GFTrackCand::addHit(unsigned int detId, unsigned int hitId, double rho, unsigned int planeId)
{
	fDetId.push_back(detId);
	fHitId.push_back(hitId);
	fPlaneId.push_back(planeId);
	fRho.push_back(rho);
}

std::vector<unsigned int> 
GFTrackCand::getHitIDs(int detId) const {
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

std::vector<unsigned int>
GFTrackCand::GetHitIDs(int detId) const {
	std::cerr << "the method GFTrackCand::GetHitIDs is deprecated. Use GFTrackCand::getHitIDs instead\n";
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

bool GFTrackCand::hitInTrack(unsigned int detId, unsigned int hitId) const
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

void GFTrackCand::Print(const Option_t* option) const {
	std::cout << "======== GFTrackCand::print ========\n";
	std::cout << "mcTrackId=" << fMcTrackId << "\n";
	std::cout << "seed values for 6D state and cov: " << std::endl;
	fState6D.Print(option);
	fCov6D.Print(option);
	std::cout << "q" << fQ << "\n";
	std::cout << "PDG code= " << fPdg << "\n";
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

void GFTrackCand::setComplTrackSeed(const TVector3& pos, const TVector3& mom, const int pdgCode, TVector3 posError, TVector3 momError){
	std::cerr << "the method GFTrackCand::setComplTrackSeed is deprecated. Use GFTrackCand::set6DSeed() or  instead\n";
	setPdgCode(pdgCode); //also sets charge
	fState6D(0,0) = pos[0];
	fState6D(1,0) = pos[1];
	fState6D(2,0) = pos[2];
	fState6D(3,0) = mom[0];
	fState6D(4,0) = mom[1];
	fState6D(5,0) = mom[2];
	fCov6D(0,0) = posError[0]*posError[0];
	fCov6D(1,1) = posError[1]*posError[1];
	fCov6D(2,2) = posError[2]*posError[2];
	fCov6D(3,3) = momError[0]*momError[0];
	fCov6D(4,4) = momError[1]*momError[1];
	fCov6D(5,5) = momError[2]*momError[2];
}


void GFTrackCand::sortHits(){
	const unsigned int nHits = getNHits(); // all 4 private vectors must have the same size.

	//a vector that will be sorted to give after sort indices to sort the other vectors
	std::vector<std::pair<double, int> > order(nHits);
	for (unsigned int i = 0; i != nHits; ++i){
		order[i] = std::make_pair(fRho[i],i);
	}
	std::stable_sort(order.begin(), order.end()); // by default sort uses the ".first" value of the pair when sorting a std container of pairs

	std::vector<unsigned int> indices(nHits);
	for (unsigned int i = 0; i != nHits; ++i){
		indices[i] = order[i].second;
	}

	sortHits(indices);
}


void GFTrackCand::sortHits(std::vector<unsigned int> indices){

	const unsigned int nHits(getNHits());
	if (indices.size() != nHits){
		abort();
		GFException exc("GFTrackCand::sortHits ==> Size of indices != number of hits!",__LINE__,__FILE__);
		throw exc;
	}

	//these containers will hold the sorted results. They are created to avoid probably slower in-place sorting
	std::vector<unsigned int> sortedDetId(nHits);
	std::vector<unsigned int> sortedHitId(nHits);
	std::vector<unsigned int> sortedPlaneId(nHits);
	std::vector<double> sortedRho(nHits);
	for (unsigned int i = 0; i != nHits; ++i){
		unsigned int sortIndex = indices[i];
		sortedDetId[i] = fDetId[sortIndex];
		sortedHitId[i] = fHitId[sortIndex];
		sortedPlaneId[i] = fPlaneId[sortIndex];
		sortedRho[i] = fRho[sortIndex];
	}
	//write the changes back to the private data members:
	fDetId = sortedDetId;
	fHitId = sortedHitId;
	fPlaneId = sortedPlaneId;
	fRho = sortedRho;
}
