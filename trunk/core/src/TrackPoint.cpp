/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "Track.h"

#include "TrackPoint.h"
#include "Exception.h"

#include <iostream>


namespace genfit {

TrackPoint::TrackPoint() :
  sortingParameter_(0), track_(NULL), thinScatterer_(NULL)
{
  ;
}

TrackPoint::TrackPoint(Track* track) :
  sortingParameter_(0), track_(track), thinScatterer_(NULL)
{
  ;
}

TrackPoint::TrackPoint(const std::vector< genfit::AbsMeasurement* >& rawMeasurements, Track* track) :
  sortingParameter_(0), track_(track), thinScatterer_(NULL)
{
  rawMeasurements_.reserve(rawMeasurements.size());

  for (std::vector<AbsMeasurement*>::const_iterator m = rawMeasurements.begin(); m != rawMeasurements.end(); ++m) {
    addRawMeasurement(*m);
  }
}

TrackPoint::TrackPoint(AbsMeasurement* rawMeasurement, Track* track) :
  sortingParameter_(0), track_(track), thinScatterer_(NULL)
{
  addRawMeasurement(rawMeasurement);
}


TrackPoint::TrackPoint(const TrackPoint& rhs) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_), thinScatterer_(NULL)
{
  // clone rawMeasurements
  for (std::vector<AbsMeasurement*>::const_iterator it = rhs.rawMeasurements_.begin(); it != rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = (*it)->clone();
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    AbsFitterInfo* fi = it->second->clone();
    fi->setTrackPoint(this);
    setFitterInfo(fi);
  }

  if (rhs.thinScatterer_ != NULL)
    thinScatterer_.reset(new ThinScatterer(*(rhs.thinScatterer_)));
}

TrackPoint::TrackPoint(const TrackPoint& rhs,
    const std::map<const AbsTrackRep*, AbsTrackRep*>& map,
    const std::vector<const genfit::AbsTrackRep*> * repsToIgnore) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_), thinScatterer_(NULL)
{
  // clone rawMeasurements
  for (std::vector<AbsMeasurement*>::const_iterator it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* m = (*it)->clone();
    addRawMeasurement(m);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    if (repsToIgnore != NULL) {
      if (std::find(repsToIgnore->begin(), repsToIgnore->end(), it->first) != repsToIgnore->end())
        continue;
    }
    AbsFitterInfo* fi = it->second->clone();
    fi->setRep(map.at(it->first));
    fi->setTrackPoint(this);
    setFitterInfo(fi);
  }

  if (rhs.thinScatterer_ != NULL)
    thinScatterer_.reset(new ThinScatterer(*(rhs.thinScatterer_)));
}


TrackPoint& TrackPoint::operator=(TrackPoint rhs) {
  swap(rhs);

  for (std::vector<AbsMeasurement*>::const_iterator it = rawMeasurements_.begin(); it!=rawMeasurements_.end(); ++it) {
    (*it)->setTrackPoint(this);
  }

  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    it->second->setTrackPoint(this);
  }

  return *this;
}


void TrackPoint::swap(TrackPoint& other) {
  std::swap(this->sortingParameter_, other.sortingParameter_);
  std::swap(this->track_, other.track_);
  std::swap(this->rawMeasurements_, other.rawMeasurements_);
  std::swap(this->fitterInfos_, other.fitterInfos_);
  this->thinScatterer_.swap(other.thinScatterer_);
}


TrackPoint::~TrackPoint() {
  // FIXME: We definitely need some smart containers or smart pointers that
  // take care of this, but so far we haven't found a convincing
  // option (2013-07-05).
  
  for (size_t i = 0; i < rawMeasurements_.size(); ++i)
    delete rawMeasurements_[i];

  std::map< const AbsTrackRep*, AbsFitterInfo* >::iterator it;
  for (it = fitterInfos_.begin(); it != fitterInfos_.end(); ++it)
    delete it->second;
}


AbsMeasurement* TrackPoint::getRawMeasurement(int i) const {
  if (i < 0)
    i += rawMeasurements_.size();

  return rawMeasurements_.at(i);
}


std::vector< AbsFitterInfo* > TrackPoint::getFitterInfos() const {
  std::vector< AbsFitterInfo* > retVal;

  if (fitterInfos_.empty())
    return retVal;

  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    retVal.push_back(it->second);
  }

  return retVal;
}


AbsFitterInfo* TrackPoint::getFitterInfo(const AbsTrackRep* rep) const {
  if (rep == NULL) {
    return fitterInfos_.at(track_->getCardinalRep());
  }
  return fitterInfos_.at(rep);
}



void TrackPoint::deleteRawMeasurements() {
  for (size_t i = 0; i < rawMeasurements_.size(); ++i)
    delete rawMeasurements_[i];

  rawMeasurements_.clear();
}


void TrackPoint::setFitterInfo(genfit::AbsFitterInfo* fitterInfo) {
  assert (fitterInfo != NULL);
  if (hasFitterInfo(fitterInfo->getRep()))
    delete fitterInfos_[fitterInfo->getRep()];

  fitterInfos_[fitterInfo->getRep()] = fitterInfo;
}


void TrackPoint::Print(const Option_t*) const {
  std::cout << "genfit::TrackPoint, belonging to Track " << track_ << "; sorting parameter = " << sortingParameter_ << "\n";
  std::cout << "contains " << rawMeasurements_.size() << " rawMeasurements and " << getFitterInfos().size() << " fitterInfos for " << fitterInfos_.size() << " TrackReps.\n";

  for (unsigned int i=0; i<rawMeasurements_.size(); ++i) {
    std::cout << "RawMeasurement Nr. " << i << "\n";
    rawMeasurements_[i]->Print();
    std::cout << "............\n";
  }

  for (std::map< const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    std::cout << "FitterInfo for TrackRep " << it->first << "\n";
    it->second->Print();
    std::cout << "............\n";
  }

  if (thinScatterer_)
    thinScatterer_->Print();

}


//
// This is modified from the auto-generated Streamer.
//
void TrackPoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::TrackPoint.
   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::TrackPoint thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //TObject::Streamer(R__b);
      R__b >> sortingParameter_;
      {
        std::vector<genfit::AbsMeasurement*,std::allocator<genfit::AbsMeasurement*> > &R__stl =  rawMeasurements_;
        R__stl.clear();
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsMeasurement));
        if (R__tcl1==0) {
          Error("rawMeasurements_ streamer","Missing the TClass object for genfit::AbsMeasurement!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        R__stl.reserve(R__n);
        for (R__i = 0; R__i < R__n; R__i++) {
          genfit::AbsMeasurement* R__t = 0;
          R__b >> R__t;
          R__stl.push_back(R__t);
        }
      }
      track_ = NULL;
      size_t nTrackReps;
      R__b >> nTrackReps;
      vFitterInfos_.resize(nTrackReps);
      for (size_t i = 0; i < nTrackReps; ++i)  {
        int id;
        R__b >> id;
        AbsFitterInfo* p = 0;
        R__b >> p;
        vFitterInfos_[id] = p;
      }
      thinScatterer_.reset();
      char flag;
      R__b >> flag;
      if (flag) {
        genfit::ThinScatterer *scatterer = 0;
        R__b >> scatterer;
        thinScatterer_.reset(new ThinScatterer(*scatterer));
      }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());


      // Fixup ownerships.
      for (size_t i = 0; i < rawMeasurements_.size(); ++i) {
        rawMeasurements_[i]->setTrackPoint(this);
      }
      for (size_t i = 0; i < vFitterInfos_.size(); ++i) {
        // May not have FitterInfos for all reps.
        if (vFitterInfos_[i])
          vFitterInfos_[i]->setTrackPoint(this);
      }
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //TObject::Streamer(R__b);
      R__b << sortingParameter_;
      {
        std::vector<genfit::AbsMeasurement*,std::allocator<genfit::AbsMeasurement*> > &R__stl =  rawMeasurements_;
        int R__n=(&R__stl) ? int(R__stl.size()) : 0;
        R__b << R__n;
        if(R__n) {
          std::vector<genfit::AbsMeasurement*,std::allocator<genfit::AbsMeasurement*> >::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
          }
        }
      }
      R__b << fitterInfos_.size();
      for (std::map<const AbsTrackRep*, AbsFitterInfo*>::const_iterator it = fitterInfos_.begin();
          it != fitterInfos_.end(); ++it)
      {
        int id = track_->getIdForRep(it->first);
        R__b << id;
        R__b << it->second;
      }
      if (thinScatterer_) {
        R__b << (char)1;
        R__b << thinScatterer_.get();
      } else {
        R__b << (char)0;
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}


void TrackPoint::fixupRepsForReading()
{
  for (size_t i = 0; i < vFitterInfos_.size(); ++i) {
    // The vector is filled such that i corresponds to the id of the TrackRep.
    
    // May not have FitterInfos for all reps.
    if (!vFitterInfos_[i])
      continue;
    fitterInfos_[track_->getTrackRep(i)] = vFitterInfos_[i];
    fitterInfos_[track_->getTrackRep(i)]->setRep(track_->getTrackRep(i));
  }
  vFitterInfos_.clear();
}

} /* End of namespace genfit */
