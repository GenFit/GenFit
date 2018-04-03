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

#include "Track.h"

#include "Exception.h"
#include "IO.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"
#include "PlanarMeasurement.h"
#include "AbsMeasurement.h"

#include "WireTrackCandHit.h"

#include <algorithm>
#include <map>

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TBuffer.h>

//#include <glog/logging.h>

//#define DEBUG


namespace genfit {

Track::Track() :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(6), covSeed_(6)
{
  ;
}


Track::Track(const TrackCand& trackCand, const MeasurementFactory<AbsMeasurement>& factory, AbsTrackRep* rep) :
  cardinalRep_(0), fitStatuses_(), stateSeed_(6), covSeed_(6)
{

  if (rep != nullptr)
    addTrackRep(rep);

  createMeasurements(trackCand, factory);

  // Copy seed information from candidate
  timeSeed_ = trackCand.getTimeSeed();
  stateSeed_ = trackCand.getStateSeed();
  covSeed_ = trackCand.getCovSeed();

  mcTrackId_ = trackCand.getMcTrackId();

  // fill cache
  fillPointsWithMeasurement();

  checkConsistency();
}

void
Track::createMeasurements(const TrackCand& trackCand, const MeasurementFactory<AbsMeasurement>& factory)
{
  // create the measurements using the factory.
  std::vector <AbsMeasurement*> factoryHits = factory.createMany(trackCand);

  if (factoryHits.size() != trackCand.getNHits()) {
    Exception exc("Track::Track ==> factoryHits.size() != trackCand->getNHits()",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // create TrackPoints
  for (unsigned int i=0; i<factoryHits.size(); ++i){
    TrackPoint* tp = new TrackPoint(factoryHits[i], this);
    tp->setSortingParameter(trackCand.getHit(i)->getSortingParameter());
    insertPoint(tp);
  }
}


Track::Track(AbsTrackRep* trackRep, const TVectorD& stateSeed) :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(stateSeed),
  covSeed_(TMatrixDSym::kUnit, TMatrixDSym(6))
{
  addTrackRep(trackRep);
}


Track::Track(AbsTrackRep* trackRep, const TVector3& posSeed, const TVector3& momSeed) :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(6),
  covSeed_(TMatrixDSym::kUnit, TMatrixDSym(6))
{
  addTrackRep(trackRep);
  setStateSeed(posSeed, momSeed);
}


Track::Track(AbsTrackRep* trackRep, const TVectorD& stateSeed, const TMatrixDSym& covSeed) :
  cardinalRep_(0), fitStatuses_(), mcTrackId_(-1), timeSeed_(0), stateSeed_(stateSeed),
  covSeed_(covSeed)
{
  addTrackRep(trackRep);
}


Track::Track(const Track& rhs) :
  TObject(rhs),
  cardinalRep_(rhs.cardinalRep_), mcTrackId_(rhs.mcTrackId_), timeSeed_(rhs.timeSeed_),
  stateSeed_(rhs.stateSeed_), covSeed_(rhs.covSeed_)
{
  rhs.checkConsistency();

  std::map<const AbsTrackRep*, AbsTrackRep*> oldRepNewRep;

  for (std::vector<AbsTrackRep*>::const_iterator it=rhs.trackReps_.begin(); it!=rhs.trackReps_.end(); ++it) {
    AbsTrackRep* newRep = (*it)->clone();
    addTrackRep(newRep);
    oldRepNewRep[(*it)] = newRep;
  }

  trackPoints_.reserve(rhs.trackPoints_.size());
  for (std::vector<TrackPoint*>::const_iterator it=rhs.trackPoints_.begin(); it!=rhs.trackPoints_.end(); ++it) {
    trackPoints_.push_back(new TrackPoint(**it, oldRepNewRep));
    trackPoints_.back()->setTrack(this);
  }

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=rhs.fitStatuses_.begin(); it!=rhs.fitStatuses_.end(); ++it) {
    setFitStatus(it->second->clone(), oldRepNewRep[it->first]);
  }

  fillPointsWithMeasurement();

  checkConsistency();
}

Track& Track::operator=(Track other) {
  swap(other);

  for (std::vector<TrackPoint*>::const_iterator it=trackPoints_.begin(); it!=trackPoints_.end(); ++it) {
    trackPoints_.back()->setTrack(this);
  }

  fillPointsWithMeasurement();

  checkConsistency();

  return *this;
}

void Track::swap(Track& other) {
  std::swap(this->trackReps_, other.trackReps_);
  std::swap(this->cardinalRep_, other.cardinalRep_);
  std::swap(this->trackPoints_, other.trackPoints_);
  std::swap(this->trackPointsWithMeasurement_, other.trackPointsWithMeasurement_);
  std::swap(this->fitStatuses_, other.fitStatuses_);
  std::swap(this->mcTrackId_, other.mcTrackId_);
  std::swap(this->timeSeed_, other.timeSeed_);
  std::swap(this->stateSeed_, other.stateSeed_);
  std::swap(this->covSeed_, other.covSeed_);

}

Track::~Track() {
  this->Clear();
}

void Track::Clear(Option_t*)
{
  // This function is needed for TClonesArray embedding.
  // FIXME: smarter containers or pointers needed ...
  for (size_t i = 0; i < trackPoints_.size(); ++i)
    delete trackPoints_[i];

  trackPoints_.clear();
  trackPointsWithMeasurement_.clear();

  for (std::map< const AbsTrackRep*, FitStatus* >::iterator it = fitStatuses_.begin(); it!= fitStatuses_.end(); ++it)
    delete it->second;
  fitStatuses_.clear();

  for (size_t i = 0; i < trackReps_.size(); ++i)
    delete trackReps_[i];
  trackReps_.clear();

  cardinalRep_ = 0;

  mcTrackId_ = -1;

  timeSeed_ = 0;
  stateSeed_.Zero();
  covSeed_.Zero();
}


TrackPoint* Track::getPoint(int id) const {
  if (id < 0)
    id += trackPoints_.size();

  return trackPoints_.at(id);
}


TrackPoint* Track::getPointWithMeasurement(int id) const {
  if (id < 0)
    id += trackPointsWithMeasurement_.size();

  return trackPointsWithMeasurement_.at(id);
}


TrackPoint* Track::getPointWithMeasurementAndFitterInfo(int id, const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (id >= 0) {
    int i = 0;
    for (std::vector<TrackPoint*>::const_iterator it = trackPointsWithMeasurement_.begin(); it != trackPointsWithMeasurement_.end(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        ++i;
      }
    }
  } else {
    // Search backwards.
    int i = -1;
    for (std::vector<TrackPoint*>::const_reverse_iterator it = trackPointsWithMeasurement_.rbegin(); it != trackPointsWithMeasurement_.rend(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        --i;
      }
    }
  }

  // Not found, i.e. abs(id) > number of fitted TrackPoints
  return 0;
}


TrackPoint* Track::getPointWithFitterInfo(int id, const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (id >= 0) {
    int i = 0;
    for (std::vector<TrackPoint*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        ++i;
      }
    }
  } else {
    // Search backwards.
    int i = -1;
    for (std::vector<TrackPoint*>::const_reverse_iterator it = trackPoints_.rbegin(); it != trackPoints_.rend(); ++it) {
      if ((*it)->hasFitterInfo(rep)) {
        if (id == i)
          return (*it);
        --i;
      }
    }
  }

  // Not found, i.e. abs(id) > number of fitted TrackPoints
  return 0;
}


const MeasuredStateOnPlane& Track::getFittedState(int id, const AbsTrackRep* rep, bool biased) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  TrackPoint* point = getPointWithFitterInfo(id, rep);
  if (point == nullptr) {
    Exception exc("Track::getFittedState ==> no trackPoint with fitterInfo for rep",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  return point->getFitterInfo(rep)->getFittedState(biased);
}


int Track::getIdForRep(const AbsTrackRep* rep) const
{
  for (size_t i = 0; i < trackReps_.size(); ++i)
    if (trackReps_[i] == rep)
      return i;

  Exception exc("Track::getIdForRep ==> cannot find TrackRep in Track",__LINE__,__FILE__);
  exc.setFatal();
  throw exc;
}


bool Track::hasFitStatus(const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (fitStatuses_.find(rep) == fitStatuses_.end())
    return false;

  return (fitStatuses_.at(rep) != nullptr);
}


bool Track::hasKalmanFitStatus(const AbsTrackRep* rep) const {
  if (rep == nullptr)
    rep = getCardinalRep();

  if (fitStatuses_.find(rep) == fitStatuses_.end())
    return false;

  return (dynamic_cast<KalmanFitStatus*>(fitStatuses_.at(rep)) != nullptr);
}


KalmanFitStatus* Track::getKalmanFitStatus(const AbsTrackRep* rep) const {
  return dynamic_cast<KalmanFitStatus*>(getFitStatus(rep));
}


void Track::setFitStatus(FitStatus* fitStatus, const AbsTrackRep* rep) {
  if (fitStatuses_.find(rep) != fitStatuses_.end())
    delete fitStatuses_.at(rep);

  fitStatuses_[rep] = fitStatus;
}


void Track::setStateSeed(const TVector3& pos, const TVector3& mom) {
  stateSeed_.ResizeTo(6);

  stateSeed_(0) = pos.X();
  stateSeed_(1) = pos.Y();
  stateSeed_(2) = pos.Z();

  stateSeed_(3) = mom.X();
  stateSeed_(4) = mom.Y();
  stateSeed_(5) = mom.Z();
}



void Track::insertPoint(TrackPoint* point, int id) {

  point->setTrack(this);

  #ifdef DEBUG
  debugOut << "Track::insertPoint at position " << id  << "\n";
  #endif
  assert(point!=nullptr);
  trackHasChanged();

  point->setTrack(this);

  if (trackPoints_.size() == 0) {
    trackPoints_.push_back(point);

    if (point->hasRawMeasurements())
      trackPointsWithMeasurement_.push_back(point);

    return;
  }

  if (id == -1 || id == (int)trackPoints_.size()) {
    trackPoints_.push_back(point);

    if (point->hasRawMeasurements())
      trackPointsWithMeasurement_.push_back(point);

    deleteReferenceInfo(std::max(0, (int)trackPoints_.size()-2), (int)trackPoints_.size()-1);

    // delete fitter infos if inserted point has a measurement
    if (point->hasRawMeasurements()) {
      deleteForwardInfo(-1, -1);
      deleteBackwardInfo(0, -2);
    }

    return;
  }

  // [-size, size-1] is the allowed range
  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size() + 1;

  // insert
  trackPoints_.insert(trackPoints_.begin() + id, point);  // insert inserts BEFORE

  // delete fitter infos if inserted point has a measurement
  if (point->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id);
  }

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id+1));

  fillPointsWithMeasurement();
}


void Track::insertPoints(std::vector<TrackPoint*> points, int id) {

  int nBefore = getNumPoints();
  int n = points.size();

  if (n == 0)
    return;
  if (n == 1) {
    insertPoint(points[0], id);
    return;
  }

  for (std::vector<TrackPoint*>::iterator p = points.begin(); p != points.end(); ++p)
    (*p)->setTrack(this);

  if (id == -1 || id == (int)trackPoints_.size()) {
    trackPoints_.insert(trackPoints_.end(), points.begin(), points.end());

    deleteReferenceInfo(std::max(0, nBefore-1), nBefore);

    deleteForwardInfo(nBefore, -1);
    deleteBackwardInfo(0, std::max(0, nBefore-1));

    fillPointsWithMeasurement();

    return;
  }


  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size() + 1;


  // insert
  trackPoints_.insert(trackPoints_.begin() + id, points.begin(), points.end());  // insert inserts BEFORE

  // delete fitter infos if inserted point has a measurement
  deleteForwardInfo(id, -1);
  deleteBackwardInfo(0, id+n);

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id));
  deleteReferenceInfo(std::max(0, id+n-1), std::min((int)trackPoints_.size()-1, id+n));

  fillPointsWithMeasurement();
}


void Track::deletePoint(int id) {

  #ifdef DEBUG
  debugOut << "Track::deletePoint at position " << id  << "\n";
  #endif

  trackHasChanged();

  if (id < 0)
    id += trackPoints_.size();
  assert(id>0);


  // delete forwardInfo after point (backwardInfo before point) if deleted point has a measurement
  if (trackPoints_[id]->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id-1);
  }

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id+1));

  // delete point
  std::vector<TrackPoint*>::iterator it = std::find(trackPointsWithMeasurement_.begin(), trackPointsWithMeasurement_.end(), trackPoints_[id]);
  if (it != trackPointsWithMeasurement_.end())
    trackPointsWithMeasurement_.erase(it);

  delete trackPoints_[id];
  trackPoints_.erase(trackPoints_.begin()+id);

  fillPointsWithMeasurement();

}


void Track::insertMeasurement(AbsMeasurement* measurement, int id) {
  insertPoint(new TrackPoint(measurement, this), id);
}
  
void Track::deleteFittedState(const genfit::AbsTrackRep* rep) {
  if(hasFitStatus(rep)) {
    delete fitStatuses_.at(rep);
    fitStatuses_.erase(rep);
  }

  // delete FitterInfos related to the deleted TrackRep
  for (const auto& trackPoint : trackPoints_) {
    if(trackPoint->hasFitterInfo(rep)) {
      trackPoint->deleteFitterInfo(rep);
    }
  }
}


void Track::mergeTrack(const Track* other, int id) {

  #ifdef DEBUG
  debugOut << "Track::mergeTrack\n";
  #endif

  if (other->getNumPoints() == 0)
    return;

  std::map<const AbsTrackRep*, AbsTrackRep*> otherRepThisRep;
  std::vector<const AbsTrackRep*> otherRepsToRemove;

  for (std::vector<AbsTrackRep*>::const_iterator otherRep=other->trackReps_.begin(); otherRep!=other->trackReps_.end(); ++otherRep) {
    bool found(false);
    for (std::vector<AbsTrackRep*>::const_iterator thisRep=trackReps_.begin(); thisRep!=trackReps_.end(); ++thisRep) {
      if ((*thisRep)->isSame(*otherRep)) {
        otherRepThisRep[*otherRep] = *thisRep;
        #ifdef DEBUG
        debugOut << " map other rep " << *otherRep << " to " << (*thisRep) << "\n";
        #endif
        if (found) {
          Exception exc("Track::mergeTrack ==> more than one matching rep.",__LINE__,__FILE__);
          exc.setFatal();
          throw exc;
        }
        found = true;
        break;
      }
    }
    if (!found) {
      otherRepsToRemove.push_back(*otherRep);
      #ifdef DEBUG
      debugOut << " remove other rep " << *otherRep << "\n";
      #endif
    }
  }


  std::vector<TrackPoint*> points;
  points.reserve(other->getNumPoints());

  for (std::vector<TrackPoint*>::const_iterator otherTp=other->trackPoints_.begin(); otherTp!=other->trackPoints_.end(); ++otherTp) {
    points.push_back(new TrackPoint(**otherTp, otherRepThisRep, &otherRepsToRemove));
  }

  insertPoints(points, id);
}


void Track::addTrackRep(AbsTrackRep* trackRep) {
  trackReps_.push_back(trackRep);
  fitStatuses_[trackRep] = new FitStatus();
}


void Track::deleteTrackRep(int id) {
  if (id < 0)
    id += trackReps_.size();

  AbsTrackRep* rep = trackReps_.at(id);

  // update cardinalRep_
  if (int(cardinalRep_) == id)
    cardinalRep_ = 0; // reset
  else if (int(cardinalRep_) > id)
    --cardinalRep_; // make cardinalRep_ point to the same TrackRep before and after deletion

  // delete FitterInfos related to the deleted TrackRep
  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin(); pointIt != trackPoints_.end(); ++pointIt) {
    (*pointIt)->deleteFitterInfo(rep);
  }

  // delete fitStatus
  delete fitStatuses_.at(rep);
  fitStatuses_.erase(rep);

  // delete rep
  delete rep;
  trackReps_.erase(trackReps_.begin()+id);
}


void Track::setCardinalRep(int id) {

  if (id < 0)
    id += trackReps_.size();

  if (id >= 0 && (unsigned int)id < trackReps_.size())
    cardinalRep_ = id;
  else {
    cardinalRep_ = 0;
    errorOut << "Track::setCardinalRep: Attempted to set cardinalRep_ to a value out of bounds. Resetting  cardinalRep_ to 0." << std::endl;
  }
}


void Track::determineCardinalRep() {

  // Todo: test

  if (trackReps_.size() <= 1)
    return;

  double minChi2(9.E99);
  const AbsTrackRep* bestRep(nullptr);

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it = fitStatuses_.begin(); it != fitStatuses_.end(); ++it) {
    if (it->second->isFitConverged()) {
      if (it->second->getChi2() < minChi2) {
        minChi2 = it->second->getChi2();
        bestRep = it->first;
      }
    }
  }

  if (bestRep != nullptr) {
    setCardinalRep(getIdForRep(bestRep));
  }
}


bool Track::sort() {
  #ifdef DEBUG
  debugOut << "Track::sort \n";
  #endif

  int nPoints(trackPoints_.size());
  // original order
  std::vector<TrackPoint*> pointsBefore(trackPoints_);

  // sort
  std::stable_sort(trackPoints_.begin(), trackPoints_.end(), TrackPointComparator());

  // see where order changed
  int equalUntil(-1), equalFrom(nPoints);
  for (int i = 0; i<nPoints; ++i) {
    if (pointsBefore[i] == trackPoints_[i])
      equalUntil = i;
    else
      break;
  }

  if (equalUntil == nPoints-1)
    return false; // sorting did not change anything


  trackHasChanged();

  for (int i = nPoints-1; i>equalUntil; --i) {
    if (pointsBefore[i] == trackPoints_[i])
      equalFrom = i;
    else
      break;
  }

  #ifdef DEBUG
  debugOut << "Track::sort. Equal up to (including) hit " << equalUntil << " and from (including) hit " << equalFrom << " \n";
  #endif

  deleteForwardInfo(equalUntil+1, -1);
  deleteBackwardInfo(0, equalFrom-1);

  deleteReferenceInfo(std::max(0, equalUntil+1), std::min((int)trackPoints_.size()-1, equalFrom-1));

  fillPointsWithMeasurement();

  return true;
}


bool Track::udpateSeed(int id, AbsTrackRep* rep, bool biased) {
  try {
    const MeasuredStateOnPlane& fittedState = getFittedState(id, rep, biased);
    setTimeSeed(fittedState.getTime());
    setStateSeed(fittedState.get6DState());
    setCovSeed(fittedState.get6DCov());

    double fittedCharge = fittedState.getCharge();

    for (unsigned int i = 0; i<trackReps_.size(); ++i) {
      if (trackReps_[i]->getPDGCharge() * fittedCharge < 0) {
        trackReps_[i]->switchPDGSign();
      }
    }
  }
  catch (Exception& e) {
    // in this case the original track seed will be used
    return false;
  }
  return true;
}


void Track::reverseTrackPoints() {

  std::reverse(trackPoints_.begin(),trackPoints_.end());

  deleteForwardInfo(0, -1);
  deleteBackwardInfo(0, -1);
  deleteReferenceInfo(0, -1);

  fillPointsWithMeasurement();
}


void Track::switchPDGSigns(AbsTrackRep* rep) {
  if (rep != nullptr) {
    rep->switchPDGSign();
    return;
  }

  for (unsigned int i = 0; i<trackReps_.size(); ++i) {
    trackReps_[i]->switchPDGSign();
  }
}


void Track::reverseTrack() {
  udpateSeed(-1); // set fitted state of last hit as new seed
  reverseMomSeed(); // flip momentum direction
  switchPDGSigns();
  reverseTrackPoints(); // also deletes all fitterInfos
}


void Track::deleteForwardInfo(int startId, int endId, const AbsTrackRep* rep) {
  #ifdef DEBUG
  debugOut << "Track::deleteForwardInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteForwardInfo();
    }
    else {
      const std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteForwardInfo();
      }
    }
  }
}

void Track::deleteBackwardInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteBackwardInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);


  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteBackwardInfo();
    }
    else {
      const std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteBackwardInfo();
      }
    }
  }
}

void Track::deleteReferenceInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteReferenceInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteReferenceInfo();
    }
    else {
      std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteReferenceInfo();
      }
    }
  }
}

void Track::deleteMeasurementInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteMeasurementInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteMeasurementInfo();
    }
    else {
      std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteMeasurementInfo();
      }
    }
  }
}

void Track::deleteFitterInfo(int startId, int endId, const AbsTrackRep* rep) {

  #ifdef DEBUG
  debugOut << "Track::deleteFitterInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != nullptr) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->deleteFitterInfo(rep);
    }
    else {
      for (std::vector<AbsTrackRep*>::const_iterator repIt = trackReps_.begin(); repIt != trackReps_.end(); ++repIt) {
        if ((*pointIt)->hasFitterInfo(*repIt))
          (*pointIt)->deleteFitterInfo(*repIt);
      }
    }
  }
}


double Track::getTrackLen(AbsTrackRep* rep, int startId, int endId) const {

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();

  bool backwards(false);
  if (startId > endId) {
    double temp = startId;
    startId = endId;
    endId = temp;
    backwards = true;
  }

  endId += 1;

  if (rep == nullptr)
    rep = getCardinalRep();

  double trackLen(0);
  StateOnPlane state;

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (! (*pointIt)->hasFitterInfo(rep)) {
      Exception e("Track::getTracklength: trackPoint has no fitterInfo", __LINE__,__FILE__);
      throw e;
    }

    if (pointIt != trackPoints_.begin() + startId) {
      trackLen += rep->extrapolateToPlane(state, (*pointIt)->getFitterInfo(rep)->getPlane());
    }

    state = (*pointIt)->getFitterInfo(rep)->getFittedState();
  }

  if (backwards)
    trackLen *= -1.;

  return trackLen;
}


TrackCand* Track::constructTrackCand() const {
  TrackCand* cand = new TrackCand();

  cand->setTime6DSeedAndPdgCode(timeSeed_, stateSeed_, getCardinalRep()->getPDG());
  cand->setCovSeed(covSeed_);
  cand->setMcTrackId(mcTrackId_);

  for (unsigned int i = 0; i < trackPointsWithMeasurement_.size(); ++i) {
    const TrackPoint* tp = trackPointsWithMeasurement_[i];
    const std::vector< AbsMeasurement* >& measurements = tp->getRawMeasurements();

    for (unsigned int j = 0; j < measurements.size(); ++j) {
      const AbsMeasurement* m = measurements[j];
      TrackCandHit* tch;

      int planeId = -1;
      if (dynamic_cast<const PlanarMeasurement*>(m)) {
        planeId = static_cast<const PlanarMeasurement*>(m)->getPlaneId();
      }

      if (m->isLeftRightMeasurement()) {
        tch = new WireTrackCandHit(m->getDetId(), m->getHitId(), planeId,
            tp->getSortingParameter(), m->getLeftRightResolution());
      }
      else {
        tch = new TrackCandHit(m->getDetId(), m->getHitId(), planeId,
            tp->getSortingParameter());
      }
      cand->addHit(tch);
    }
  }

  return cand;
}


double Track::getTOF(AbsTrackRep* rep, int startId, int endId) const {

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();

  if (startId > endId) {
    std::swap(startId, endId);
  }

  endId += 1;

  if (rep == nullptr)
    rep = getCardinalRep();

  StateOnPlane state;

  const TrackPoint* startPoint(trackPoints_[startId]);
  const TrackPoint* endPoint(trackPoints_[endId]);
  
  if (!startPoint->hasFitterInfo(rep)
      || !endPoint->hasFitterInfo(rep)) {
      Exception e("Track::getTOF: trackPoint has no fitterInfo", __LINE__,__FILE__);
      throw e;
    }

  double tof = (rep->getTime(endPoint->getFitterInfo(rep)->getFittedState())
		- rep->getTime(startPoint->getFitterInfo(rep)->getFittedState()));
  return tof;
}


void Track::fixWeights(AbsTrackRep* rep, int startId, int endId) {

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();

  assert(startId >= 0);
  assert(startId <= endId);
  assert(endId <= (int)trackPoints_.size());

  std::vector< AbsFitterInfo* > fis;

  for (std::vector<TrackPoint*>::iterator tp = trackPoints_.begin() + startId; tp != trackPoints_.begin() + endId; ++tp) {
    fis.clear();
    if (rep == nullptr) {
      fis = (*tp)->getFitterInfos();
    }
    else if ((*tp)->hasFitterInfo(rep)) {
      fis.push_back((*tp)->getFitterInfo(rep));
    }

    for (std::vector< AbsFitterInfo* >::iterator fi = fis.begin(); fi != fis.end(); ++fi) {
      KalmanFitterInfo* kfi = dynamic_cast<KalmanFitterInfo*>(*fi);
      if (kfi == nullptr)
        continue;

      kfi->fixWeights();
    }
  }
}


void Track::prune(const Option_t* option) {

  PruneFlags f;
  f.setFlags(option);

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->getPruneFlags().setFlags(option);
  }

  // prune trackPoints
  if (f.hasFlags("F") || f.hasFlags("L")) {
    TrackPoint* firstPoint = getPointWithFitterInfo(0);
    TrackPoint* lastPoint = getPointWithFitterInfo(-1);
    for (unsigned int i = 0; i<trackPoints_.size(); ++i) {
      if (trackPoints_[i] == firstPoint && f.hasFlags("F"))
        continue;

      if (trackPoints_[i] == lastPoint && f.hasFlags("L"))
        continue;

      delete trackPoints_[i];
      trackPoints_.erase(trackPoints_.begin()+i);
      --i;
    }
  }

  // prune TrackReps
  if (f.hasFlags("C")) {
    for (unsigned int i = 0; i < trackReps_.size(); ++i) {
      if (i != cardinalRep_) {
        deleteTrackRep(i);
        --i;
      }
    }
  }


  // from remaining trackPoints: prune measurementsOnPlane, unneeded fitterInfoStuff
  for (unsigned int i = 0; i<trackPoints_.size(); ++i) {
    if (f.hasFlags("W"))
      trackPoints_[i]->deleteRawMeasurements();

    std::vector< AbsFitterInfo* > fis =  trackPoints_[i]->getFitterInfos();
    for (unsigned int j = 0; j<fis.size(); ++j) {

      if (i == 0 && f.hasFlags("FLI"))
        fis[j]->deleteForwardInfo();
      else if (i == trackPoints_.size()-1 && f.hasFlags("FLI"))
        fis[j]->deleteBackwardInfo();
      else if (f.hasFlags("FI"))
        fis[j]->deleteForwardInfo();
      else if (f.hasFlags("LI"))
        fis[j]->deleteBackwardInfo();

      if (f.hasFlags("U") && dynamic_cast<KalmanFitterInfo*>(fis[j]) != nullptr) {
        static_cast<KalmanFitterInfo*>(fis[j])->deletePredictions();
      }

      // also delete reference info if points have been removed since it is invalid then!
      if (f.hasFlags("R") or f.hasFlags("F") or f.hasFlags("L"))
        fis[j]->deleteReferenceInfo();
      if (f.hasFlags("M"))
        fis[j]->deleteMeasurementInfo();
    }
  }

  fillPointsWithMeasurement();

  #ifdef DEBUG
  debugOut << "pruned Track: "; Print();
  #endif

}


void Track::Print(const Option_t* option) const {
  TString opt = option;
  opt.ToUpper();
  if (opt.Contains("C")) { // compact

    printOut << "\n    ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {

      int color = 32*(size_t)(trackPoints_[i]) % 15;
      switch (color) {
        case 0:
          printOut<<"\033[1;30m";
          break;
        case 1:
          printOut<<"\033[0;34m";
          break;
        case 2:
          printOut<<"\033[1;34m";
          break;
        case 3:
          printOut<<"\033[0;32m";
          break;
        case 4:
          printOut<<"\033[1;32m";
          break;
        case 5:
          printOut<<"\033[0;36m";
          break;
        case 6:
          printOut<<"\033[1;36m";
          break;
        case 7:
          printOut<<"\033[0;31m";
          break;
        case 8:
          printOut<<"\033[1;31m";
          break;
        case 9:
          printOut<<"\033[0;35m";
          break;
        case 10:
          printOut<<"\033[1;35m";
          break;
        case 11:
          printOut<<"\033[0;33m";
          break;
        case 12:
          printOut<<"\033[1;33m";
          break;
        case 13:
          printOut<<"\033[0;37m";
          break;
        default:
          ;
      }
      printOut << trackPoints_[i] << "\033[00m  ";
    }
    printOut << "\n";

    printOut << "   ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {
      printf("% -9.3g  ", trackPoints_[i]->getSortingParameter());
    }

    for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
      printOut << "\n" << getIdForRep(*rep) << "   ";
      for (unsigned int i=0; i<trackPoints_.size(); ++i) {
        if (! trackPoints_[i]->hasFitterInfo(*rep)) {
          printOut << "           ";
          continue;
        }
        AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(*rep);
        if (fi->hasMeasurements())
          printOut << "M";
        else
          printOut << " ";

        if (fi->hasReferenceState())
          printOut << "R";
        else
          printOut << " ";

        printOut << "         ";
      }
      printOut << "\n";

      printOut << " -> ";
      for (unsigned int i=0; i<trackPoints_.size(); ++i) {
        if (! trackPoints_[i]->hasFitterInfo(*rep)) {
          printOut << "           ";
          continue;
        }
        AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(*rep);
        if (fi->hasForwardPrediction())
          printOut << "P";
        else
          printOut << " ";

        if (fi->hasForwardUpdate())
          printOut << "U";
        else
          printOut << " ";

        printOut << "         ";
      }
      printOut << "\n";

      printOut << " <- ";
      for (unsigned int i=0; i<trackPoints_.size(); ++i) {
        if (! trackPoints_[i]->hasFitterInfo(*rep)) {
          printOut << "           ";
          continue;
        }
        AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(*rep);
        if (fi->hasBackwardPrediction())
          printOut << "P";
        else
          printOut << " ";

        if (fi->hasBackwardUpdate())
          printOut << "U";
        else
          printOut << " ";

        printOut << "         ";
      }

      printOut << "\n";

    } //end loop over reps

    printOut << "\n";
    return;
  }



  printOut << "=======================================================================================\n";
  printOut << "genfit::Track, containing " << trackPoints_.size() << " TrackPoints and " << trackReps_.size() << " TrackReps.\n";
  printOut << " Seed state: "; stateSeed_.Print();

  for (unsigned int i=0; i<trackReps_.size(); ++i) {
    printOut << " TrackRep Nr. " << i;
    if (i == cardinalRep_)
      printOut << " (This is the cardinal rep)";
    printOut << "\n";
    trackReps_[i]->Print();
  }

  printOut << "---------------------------------------------------------------------------------------\n";

  for (unsigned int i=0; i<trackPoints_.size(); ++i) {
    printOut << "TrackPoint Nr. " << i << "\n";
    trackPoints_[i]->Print();
    printOut << "..........................................................................\n";
  }

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->Print();
  }

  printOut << "=======================================================================================\n";

}


void Track::checkConsistency() const {

  bool consistent = true;
  std::stringstream failures;

  std::map<const AbsTrackRep*, const KalmanFitterInfo*> prevFis;

  // check if seed is 6D
  if (stateSeed_.GetNrows() != 6) {
    failures << "Track::checkConsistency(): stateSeed_ dimension != 6" << std::endl;
    consistent = false;
  }

  if (covSeed_.GetNrows() != 6) {
    failures << "Track::checkConsistency(): covSeed_ dimension != 6" << std::endl;
    consistent = false;
  }

  if (covSeed_.Max() == 0.) {
    // Nota bene: The consistency is not set to false when this occurs, because it does not break the consistency of
    // the track. However, when something else fails we keep this as additional error information.
    failures << "Track::checkConsistency(): Warning: covSeed_ is zero" << std::endl;
  }

  // check if correct number of fitStatuses
  if (fitStatuses_.size() != trackReps_.size()) {
    failures << "Track::checkConsistency(): Number of fitStatuses is != number of TrackReps " << std::endl;
    consistent = false;
  }

  // check if cardinalRep_ is in range of trackReps_
  if (trackReps_.size() && cardinalRep_ >= trackReps_.size()) {
    failures << "Track::checkConsistency(): cardinalRep id " << cardinalRep_ << " out of bounds" << std::endl;
    consistent = false;
  }

  for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
    // check for nullptr
    if ((*rep) == nullptr) {
      failures << "Track::checkConsistency(): TrackRep is nullptr" << std::endl;
      consistent = false;
    }

    // check for valid pdg code
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle((*rep)->getPDG());
    if (particle == nullptr) {
      failures << "Track::checkConsistency(): TrackRep pdg ID " << (*rep)->getPDG() << " is not valid" << std::endl;
      consistent = false;
    }

    // check if corresponding FitStatus is there
    if (fitStatuses_.find(*rep) == fitStatuses_.end() and fitStatuses_.find(*rep)->second != nullptr) {
      failures << "Track::checkConsistency(): No FitStatus for Rep or FitStatus is nullptr" << std::endl;
      consistent = false;
    }
  }

  // check TrackPoints
  for (std::vector<TrackPoint*>::const_iterator tp = trackPoints_.begin(); tp != trackPoints_.end(); ++tp) {
    // check for nullptr
    if ((*tp) == nullptr) {
      failures << "Track::checkConsistency(): TrackPoint is nullptr" << std::endl;
      consistent = false;
    }
    // check if trackPoint points back to this track
    if ((*tp)->getTrack() != this) {
      failures << "Track::checkConsistency(): TrackPoint does not point back to this track" << std::endl;
      consistent = false;
    }

    // check rawMeasurements
    const std::vector<AbsMeasurement*>& rawMeasurements = (*tp)->getRawMeasurements();
    for (std::vector<AbsMeasurement*>::const_iterator m = rawMeasurements.begin(); m != rawMeasurements.end(); ++m) {
      // check for nullptr
      if ((*m) == nullptr) {
        failures << "Track::checkConsistency(): Measurement is nullptr" << std::endl;
        consistent = false;
      }
      // check if measurement points back to TrackPoint
      if ((*m)->getTrackPoint() != *tp) {
        failures << "Track::checkConsistency(): Measurement does not point back to correct TrackPoint" << std::endl;
        consistent = false;
      }
    }

    // check fitterInfos
    std::vector<AbsFitterInfo*> fitterInfos = (*tp)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fi = fitterInfos.begin(); fi != fitterInfos.end(); ++fi) {
      // check for nullptr
      if ((*fi) == nullptr) {
        failures << "Track::checkConsistency(): FitterInfo is nullptr. TrackPoint: " << *tp << std::endl;
        consistent = false;
      }

      // check if fitterInfos point to valid TrackReps in trackReps_
      int mycount (0);
      for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
        if ( (*rep) == (*fi)->getRep() ) {
          ++mycount;
        }
      }
      if (mycount ==  0) {
        failures << "Track::checkConsistency(): fitterInfo points to TrackRep which is not in Track" << std::endl;
        consistent = false;
      }

      if (!( (*fi)->checkConsistency(&(this->getFitStatus((*fi)->getRep())->getPruneFlags())) ) ) {
        failures << "Track::checkConsistency(): FitterInfo not consistent. TrackPoint: " << *tp << std::endl;
        consistent = false;
      }

      if (dynamic_cast<KalmanFitterInfo*>(*fi) != nullptr) {
        if (prevFis[(*fi)->getRep()] != nullptr &&
            static_cast<KalmanFitterInfo*>(*fi)->hasReferenceState() &&
            prevFis[(*fi)->getRep()]->hasReferenceState() ) {
          double len = static_cast<KalmanFitterInfo*>(*fi)->getReferenceState()->getForwardSegmentLength();
          double prevLen = prevFis[(*fi)->getRep()]->getReferenceState()->getBackwardSegmentLength();
          if (fabs(prevLen + len) > 1E-10 ) {
            failures << "Track::checkConsistency(): segment lengths of reference states for rep " << (*fi)->getRep() << " (id " << getIdForRep((*fi)->getRep()) << ") at TrackPoint " << (*tp) << " don't match" << std::endl;
            failures << prevLen << " + " << len << " = " << prevLen + len << std::endl;
            failures << "TrackPoint " << *tp << ", FitterInfo " << *fi << ", rep " << getIdForRep((*fi)->getRep()) << std::endl;
            consistent = false;
          }
        }

        prevFis[(*fi)->getRep()] = static_cast<KalmanFitterInfo*>(*fi);
      }
      else
        prevFis[(*fi)->getRep()] = nullptr;

    } // end loop over FitterInfos

  } // end loop over TrackPoints


  // check trackPointsWithMeasurement_
  std::vector<TrackPoint*> trackPointsWithMeasurement;

  for (std::vector<TrackPoint*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    if ((*it)->hasRawMeasurements()) {
      trackPointsWithMeasurement.push_back(*it);
    }
  }

  if (trackPointsWithMeasurement.size() != trackPointsWithMeasurement_.size()) {
    failures << "Track::checkConsistency(): trackPointsWithMeasurement_ has incorrect size" << std::endl;
    consistent = false;
  }

  for (unsigned int i = 0; i < trackPointsWithMeasurement.size(); ++i) {
    if (trackPointsWithMeasurement[i] != trackPointsWithMeasurement_[i]) {
      failures << "Track::checkConsistency(): trackPointsWithMeasurement_ is not correct" << std::endl;
      failures << "has         id " << i << ", address " << trackPointsWithMeasurement_[i] << std::endl;
      failures << "should have id " << i << ", address " << trackPointsWithMeasurement[i] << std::endl;
      consistent = false;
    }
  }

  if (not consistent) {
    throw genfit::Exception(failures.str(), __LINE__, __FILE__);
  }
}


void Track::trackHasChanged() {

  #ifdef DEBUG
  debugOut << "Track::trackHasChanged \n";
  #endif

  if (fitStatuses_.empty())
    return;

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->setHasTrackChanged();
  }
}


void Track::fillPointsWithMeasurement() {
  trackPointsWithMeasurement_.clear();
  trackPointsWithMeasurement_.reserve(trackPoints_.size());

  for (std::vector<TrackPoint*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    if ((*it)->hasRawMeasurements()) {
      trackPointsWithMeasurement_.push_back(*it);
    }
  }
}


void Track::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::Track.
  const bool streamTrackPoints = true; // debugging option
   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::Track thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      this->Clear();
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      {
        std::vector<AbsTrackRep*> &R__stl =  trackReps_;
        R__stl.clear();
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
        if (R__tcl1==0) {
          Error("trackReps_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        R__stl.reserve(R__n);
        for (R__i = 0; R__i < R__n; R__i++) {
          genfit::AbsTrackRep* R__t;
          R__b >> R__t;
          R__stl.push_back(R__t);
        }
      }
      R__b >> cardinalRep_;
      if (streamTrackPoints)
      {
        std::vector<TrackPoint*> &R__stl =  trackPoints_;
        R__stl.clear();
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::TrackPoint));
        if (R__tcl1==0) {
          Error("trackPoints_ streamer","Missing the TClass object for genfit::TrackPoint!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        R__stl.reserve(R__n);
        for (R__i = 0; R__i < R__n; R__i++) {
          genfit::TrackPoint* R__t;
          R__t = (genfit::TrackPoint*)R__b.ReadObjectAny(R__tcl1);
          R__t->setTrack(this);
          R__t->fixupRepsForReading();
          R__stl.push_back(R__t);
        }
      }
      {
        std::map<const AbsTrackRep*,FitStatus*> &R__stl =  fitStatuses_;
        R__stl.clear();
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
        if (R__tcl1==0) {
          Error("fitStatuses_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
          return;
        }
        TClass *R__tcl2 = TBuffer::GetClass(typeid(genfit::FitStatus));
        if (R__tcl2==0) {
          Error("fitStatuses_ streamer","Missing the TClass object for genfit::FitStatus!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        for (R__i = 0; R__i < R__n; R__i++) {
          int id;
          R__b >> id;
          genfit::FitStatus* R__t2;
          R__t2 = (genfit::FitStatus*)R__b.ReadObjectAny(R__tcl2);

          R__stl[this->getTrackRep(id)] = R__t2;
        }
      }
      // timeSeed_ was introduced in version 3.  If reading an earlier
      // version, default to zero.
      if (R__v >= 3) {
	R__b >> timeSeed_;
      } else {
	timeSeed_ = 0;
      }
      stateSeed_.Streamer(R__b);
      covSeed_.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());

      fillPointsWithMeasurement();
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      {
        std::vector<AbsTrackRep*> &R__stl =  trackReps_;
        int R__n=int(R__stl.size());
        R__b << R__n;
        if(R__n) {
          TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
          if (R__tcl1==0) {
            Error("trackReps_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
            return;
          }
          std::vector<AbsTrackRep*>::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << *R__k;
          }
        }
      }
      R__b << cardinalRep_;
      if (streamTrackPoints)
      {
        std::vector<TrackPoint*> &R__stl =  trackPoints_;
        int R__n=int(R__stl.size());
        R__b << R__n;
        if(R__n) {
          std::vector<TrackPoint*>::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
          }
        }
      }
      {
        std::map<const AbsTrackRep*,FitStatus*> &R__stl =  fitStatuses_;
        int R__n=int(R__stl.size());
        R__b << R__n;
        if(R__n) {
          TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
          if (R__tcl1==0) {
            Error("fitStatuses_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
            return;
          }
          TClass *R__tcl2 = TBuffer::GetClass(typeid(genfit::FitStatus));
          if (R__tcl2==0) {
            Error("fitStatuses_ streamer","Missing the TClass object for genfit::FitStatus!");
            return;
          }
          std::map<const AbsTrackRep*,FitStatus*>::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            int id = this->getIdForRep((*R__k).first);
            R__b << id;
            R__b << ((*R__k).second);
          }
        }
      }
      R__b << timeSeed_;
      stateSeed_.Streamer(R__b);
      covSeed_.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

void Track::deleteTrackPointsAndFitStatus() {
  for (size_t i = 0; i < trackPoints_.size(); ++i)
    delete trackPoints_[i];

  trackPoints_.clear();
  trackPointsWithMeasurement_.clear();

  for (std::map< const AbsTrackRep*, FitStatus* >::iterator it = fitStatuses_.begin(); it!= fitStatuses_.end(); ++it)
    delete it->second;
  fitStatuses_.clear();
}

} /* End of namespace genfit */
