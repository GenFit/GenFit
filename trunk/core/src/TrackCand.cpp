/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

#include "TrackCand.h"
#include "Exception.h"
#include "TDatabasePDG.h"

#include <algorithm>
#include <iostream>
#include <utility>

namespace genfit {


TrackCand::TrackCand() :
  mcTrackId_(-1),
  pdg_(0),
  state6D_(6),
  q_(0)
{
  ;
}

TrackCand::~TrackCand() {
  for (unsigned int i=0; i<hits_.size(); ++i) {
    delete hits_[i];
  }
  hits_.clear();
}


TrackCand::TrackCand( const TrackCand& other ) :
  mcTrackId_(other.mcTrackId_),
  pdg_(other.pdg_),
  state6D_(other.state6D_),
  q_(other.q_)
{
  // deep copy
  hits_.reserve(other.hits_.size());
  for (unsigned int i=0; i<other.hits_.size(); ++i) {
    hits_.push_back( (other.hits_[i])->clone() );
  }
}

TrackCand& TrackCand::operator=(TrackCand other) {
  swap(other);
  return *this;
}


void TrackCand::swap(TrackCand& other) {
  // by swapping the members of two classes,
  // the two classes are effectively swapped
  std::swap(this->hits_, other.hits_);
  std::swap(this->mcTrackId_, other.mcTrackId_);
  std::swap(this->pdg_, other.pdg_);
  std::swap(this->state6D_, other.state6D_);
  std::swap(this->q_, other.q_);
}


TrackCandHit* TrackCand::getHit(int i) const {
  if (i < 0)
    i += hits_.size();

  return hits_.at(i);
}


void TrackCand::getHit(int i, int& detId, int& hitId) const {
  if (i < 0)
    i += hits_.size();

  detId = hits_.at(i)->getDetId();
  hitId = hits_[i]->getHitId();
}


void TrackCand::getHit(int i, int& detId, int& hitId, double& sortingParameter) const {
  if (i < 0)
    i += hits_.size();

  detId = hits_.at(i)->getDetId();
  hitId = hits_[i]->getHitId();
  sortingParameter = hits_[i]->getSortingParameter();
}


void TrackCand::getHitWithPlane(int i, int& detId, int& hitId, int& planeId) const {
  if (i < 0)
    i += hits_.size();

  detId = hits_.at(i)->getDetId();
  hitId = hits_[i]->getHitId();
  planeId = hits_[i]->getPlaneId();
}


void TrackCand::addHit(int detId, int hitId, int planeId, double sortingParameter)
{
  hits_.push_back(new TrackCandHit(detId, hitId, planeId, sortingParameter));
}

std::vector<int> TrackCand::getHitIDs(int detId) const {
  std::vector<int> result;
  for(unsigned int i=0; i<hits_.size(); ++i){
    if(detId==-2 || hits_[i]->getDetId() == detId) {
      result.push_back(hits_[i]->getHitId());
    }
  }
  return result;
}

std::vector<int> TrackCand::getDetIDs() const {
  std::vector<int> result;
  for(unsigned int i=0; i<hits_.size(); ++i){
    result.push_back(hits_[i]->getDetId());
  }
  return result;
}

std::vector<double> TrackCand::getSortingParameters() const {
  std::vector<double> result;
  for(unsigned int i=0; i<hits_.size(); ++i){
    result.push_back(hits_[i]->getSortingParameter());
  }
  return result;
}

std::set<int> TrackCand::getUniqueDetIDs() const {
  std::set<int> retVal;
  for (unsigned int i = 0; i < hits_.size(); ++i) {
    retVal.insert(hits_[i]->getDetId());
  }
  return retVal;
}


void TrackCand::setPdgCode(int pdgCode) {
  pdg_ = pdgCode;
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(pdg_);
  q_ = part->Charge() / (3.);
}


void TrackCand::reset()
{
  for (unsigned int i=0; i<hits_.size(); ++i) {
    delete hits_[i];
  }
  hits_.clear();
}


bool TrackCand::hitInTrack(int detId, int hitId) const
{
  for (unsigned int i = 0; i < hits_.size(); ++i){
    if (detId == hits_[i]->getDetId() && hitId == hits_[i]->getHitId())
      return true;
  }
  return false;
}


bool operator== (const TrackCand& lhs, const TrackCand& rhs){
  if(lhs.getNHits() != rhs.getNHits()) return false;
  for (unsigned int i = 0; i < lhs.getNHits(); ++i){
    if (lhs.getHit(i) != rhs.getHit(i)) return false;
  }
  return true;
}


void TrackCand::Print(const Option_t* option) const {
  std::cout << "======== TrackCand::print ========\n";
  std::cout << "mcTrackId=" << mcTrackId_ << "\n";
  std::cout << "seed values for 6D state: \n";
  state6D_.Print(option);
  std::cout << "q" << q_ << "\n";
  std::cout << "PDG code= " << pdg_ << "\n";
  for(unsigned int i=0; i<hits_.size(); ++i){
    hits_[i]->Print();
  }
}


void TrackCand::append(const TrackCand& rhs){
  for(unsigned int i=0; i<rhs.getNHits(); ++i){
    addHit(rhs.getHit(i)->clone());
  }
}


void TrackCand::sortHits(){
  std::stable_sort(hits_.begin(), hits_.end(), compareTrackCandHits);
}


void TrackCand::sortHits(const std::vector<unsigned int>& indices){

  const unsigned int nHits(getNHits());
  if (indices.size() != nHits){
    abort();
    Exception exc("TrackCand::sortHits ==> Size of indices != number of hits!",__LINE__,__FILE__);
    throw exc;
  }

  //these containers will hold the sorted results. They are created to avoid probably slower in-place sorting
  std::vector<TrackCandHit*> sortedHits(nHits);
  for (unsigned int i=0; i<nHits; ++i){
    sortedHits[i] = hits_[indices[i]];
  }
  //write the changes back to the private data members:
  hits_ = sortedHits;
}


void TrackCand::set6DSeed(const TVectorD& state6D, const double charge) {
  q_ = charge;
  state6D_ = state6D;
}

void TrackCand::set6DSeedAndPdgCode(const TVectorD& state6D, const int pdgCode) {
  setPdgCode(pdgCode);
  state6D_ = state6D;
}

void TrackCand::setPosMomSeed(const TVector3& pos, const TVector3& mom, const double charge) {
  q_ = charge;
  state6D_[0] = pos[0];  state6D_[1] = pos[1];  state6D_[2] = pos[2];
  state6D_[3] = mom[0];  state6D_[4] = mom[1];  state6D_[5] = mom[2];
}

void TrackCand::setPosMomSeedAndPdgCode(const TVector3& pos, const TVector3& mom, const int pdgCode) {
  setPdgCode(pdgCode);
  state6D_[0] = pos[0];  state6D_[1] = pos[1];  state6D_[2] = pos[2];
  state6D_[3] = mom[0];  state6D_[4] = mom[1];  state6D_[5] = mom[2];
}


} /* End of namespace genfit */
