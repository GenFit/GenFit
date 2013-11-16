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


#include "GFRaveConverters.h"

#include "Exception.h"

#include "rave/Plane.h"

#include "GFRaveTrackParameters.h"

#include <iostream>


namespace genfit {


std::vector < rave::Track >
GFTracksToTracks(const std::vector < genfit::Track* >  & GFTracks,
    std::map<int, genfit::trackAndState>& IdGFTrackStateMap,
    int startID){

  unsigned int ntracks(GFTracks.size());

  std::vector < rave::Track > ravetracks;
  ravetracks.reserve(ntracks);

  for (unsigned int i=0; i<ntracks; ++i){

    if (GFTracks[i] == NULL) {
      Exception exc("GFTracksToTracks ==> genfit::Track is NULL",__LINE__,__FILE__);
      throw exc;
    }

    // only convert successfully fitted tracks!
    if (!GFTracks[i]->getFitStatus(GFTracks[i]->getCardinalRep())->isFitConverged()) continue;

    if (IdGFTrackStateMap.count(startID) > 0){
      Exception exc("GFTracksToTracks ==> IdGFTrackStateMap has already an entry for this id",__LINE__,__FILE__);
      throw exc;
    }
    IdGFTrackStateMap[startID].track_ = GFTracks[i];
    IdGFTrackStateMap[startID].state_ = new MeasuredStateOnPlane(GFTracks[i]->getFittedState()); // here clones are made so that the state of the original GFTracks and their TrackReps will not be altered by the vertexing process

    ravetracks.push_back(GFTrackToTrack(IdGFTrackStateMap[startID], startID) );

    ++startID;
  }

  //std::cout << "IdGFTrackStateMap size " << IdGFTrackStateMap.size() << std::endl;
  return ravetracks;
}


rave::Track
GFTrackToTrack(trackAndState trackAndState, int id, std::string tag){

  if (trackAndState.track_ == NULL) {
    Exception exc("GFTrackToTrack ==> originaltrack is NULL",__LINE__,__FILE__);
    throw exc;
  }

  if (! trackAndState.track_->getFitStatus()->isFitConverged()) {
    Exception exc("GFTrackToTrack ==> Trackfit is not converged",__LINE__,__FILE__);
    throw exc;
  }

  TVector3 pos, mom;
  TMatrixDSym cov;

  trackAndState.track_->getFittedState().getPosMomCov(pos, mom, cov);

  // state
  rave::Vector6D ravestate(pos.X(), pos.Y(), pos.Z(),
                           mom.X(), mom.Y(), mom.Z());

  // covariance
  rave::Covariance6D ravecov(cov(0,0), cov(1,0), cov(2,0),
                             cov(1,1), cov(2,1), cov(2,2),
                             cov(3,0), cov(4,0), cov(5,0),
                             cov(3,1), cov(4,1), cov(5,1),
                             cov(3,2), cov(4,2), cov(5,2),
                             cov(3,3), cov(4,3), cov(5,3),
                             cov(4,4), cov(5,4), cov(5,5));

  //std::cerr<<"create rave track with id " << id << std::endl;
  //std::cerr<<"  pos: "; Point3DToTVector3(ravestate.position()).Print();
  //std::cerr<<"  mom: "; Vector3DToTVector3(ravestate.momentum()).Print();

  rave::Track ret(id, ravestate, ravecov,
      trackAndState.track_->getFitStatus()->getCharge(),
      trackAndState.track_->getFitStatus()->getChi2(),
      trackAndState.track_->getFitStatus()->getNdf(),
      static_cast<void*>(const_cast<Track*>(trackAndState.track_)), tag);

  //std::cout << "ret.originalObject() " << ret.originalObject() << "\n";
  //std::cout << "ret.id() " << ret.id() << "\n";

  return ret;
}


void
setData(const rave::Track& orig, MeasuredStateOnPlane* state){

  state->setPosMomCov(TVector3(orig.state().x(), orig.state().y(), orig.state().z()),
                      TVector3(orig.state().px(), orig.state().py(), orig.state().pz()),
                      Covariance6DToTMatrixDSym(orig.error()));

}


GFRaveVertex*
RaveToGFVertex(const rave::Vertex & raveVertex,
    const std::map<int, genfit::trackAndState>& IdGFTrackStateMap){

  if (!(raveVertex.isValid())) {
    Exception exc("RaveToGFVertex ==> rave Vertex is not valid!",__LINE__,__FILE__);
    throw exc;
  }

  std::vector < std::pair < float, rave::Track > > raveWeightedTracks(raveVertex.weightedTracks());
  std::vector < std::pair < float, rave::Track > > raveSmoothedTracks(raveVertex.weightedRefittedTracks());

  int id;
  unsigned int nTrks(raveWeightedTracks.size());

  // check if rave vertex has  refitted tracks
  bool smoothing(true);
  if (! (raveVertex.hasRefittedTracks()) ) {
    smoothing = false;
  }

  // check numbers of tracks and smoothed tracks
  if (smoothing && nTrks != raveSmoothedTracks.size()){
    Exception exc("RaveToGFVertex ==> number of smoothed tracks != number of tracks",__LINE__,__FILE__);
    throw exc;
  }

  // (smoothed) track parameters
  std::vector < GFRaveTrackParameters* > trackParameters;
  trackParameters.reserve(nTrks);

  // convert tracks
  for (unsigned int i=0; i<nTrks; ++i){
    id = raveWeightedTracks[i].second.id();

    if (IdGFTrackStateMap.count(id) == 0){
      Exception exc("RaveToGFVertex ==> rave track id is not present in IdGFTrackStateMap",__LINE__,__FILE__);
      throw exc;
    }

    GFRaveTrackParameters* trackparams;

    if(smoothing) {
      // convert smoothed track parameters
      trackparams = new GFRaveTrackParameters(IdGFTrackStateMap.at(id).track_, //track
          IdGFTrackStateMap.at(id).state_, //state
          raveWeightedTracks[i].first, //weight
          Vector6DToTVectorD(raveSmoothedTracks[i].second.state()), //smoothed state
          Covariance6DToTMatrixDSym(raveSmoothedTracks[i].second.error()), //smoothed cov
          true);
    }
    else {
      // convert track parameters, no smoothed tracks available
      trackparams = new GFRaveTrackParameters(IdGFTrackStateMap.at(id).track_, //track
          IdGFTrackStateMap.at(id).state_, //state
          raveWeightedTracks[i].first, //weight
          Vector6DToTVectorD(raveWeightedTracks[i].second.state()), //state
          Covariance6DToTMatrixDSym(raveWeightedTracks[i].second.error()), //cov
          false);
    }
    trackParameters.push_back(trackparams);
  }

  return new GFRaveVertex(Point3DToTVector3(raveVertex.position()),
                          Covariance3DToTMatrixDSym(raveVertex.error()),
                          trackParameters,
                          raveVertex.ndf(), raveVertex.chiSquared(), raveVertex.id());
}

void
RaveToGFVertices(std::vector<GFRaveVertex*> * GFVertices,
    const std::vector<rave::Vertex> & raveVertices,
    const std::map<int, genfit::trackAndState>& IdGFTrackStateMap){

  unsigned int nVert(raveVertices.size());

  GFVertices->reserve(nVert);

  for (unsigned int i=0; i<nVert; ++i){
    GFVertices->push_back(RaveToGFVertex(raveVertices[i], IdGFTrackStateMap));
  }
}


SharedPlanePtr
PlaneToGFDetPlane(const ravesurf::Plane& rplane) {
  return SharedPlanePtr(new DetPlane(Point3DToTVector3(rplane.position()),
                    Vector3DToTVector3(rplane.normalVector()) ));
}


TVector3
Point3DToTVector3(const rave::Point3D& v) {
  return TVector3(v.x(), v.y(), v.z());
}

TVector3
Vector3DToTVector3(const rave::Vector3D& v) {
  return TVector3(v.x(), v.y(), v.z());
}


TMatrixDSym
Covariance3DToTMatrixDSym(const rave::Covariance3D& ravecov){
  TMatrixDSym cov(3);

  cov(0,0) = ravecov.dxx();
  cov(0,1) = ravecov.dxy();
  cov(0,2) = ravecov.dxz();

  cov(1,0) = ravecov.dxy();
  cov(1,1) = ravecov.dyy();
  cov(1,2) = ravecov.dyz();

  cov(2,0) = ravecov.dxz();
  cov(2,1) = ravecov.dyz();
  cov(2,2) = ravecov.dzz();

  return cov;
}


TVectorD
Vector6DToTVectorD(const rave::Vector6D& ravevec){
  TVectorD vec(6);

  vec[0] = ravevec.x();
  vec[1] = ravevec.y();
  vec[2] = ravevec.z();

  vec[3] = ravevec.px();
  vec[4] = ravevec.py();
  vec[5] = ravevec.pz();

  return vec;
}


TMatrixDSym
Covariance6DToTMatrixDSym(const rave::Covariance6D& ravecov){
  TMatrixDSym cov(6);

  cov(0,0) = ravecov.dxx();
  cov(0,1) = ravecov.dxy();
  cov(0,2) = ravecov.dxz();
  cov(0,3) = ravecov.dxpx();
  cov(0,4) = ravecov.dxpy();
  cov(0,5) = ravecov.dxpz();

  cov(1,0) = ravecov.dxy();
  cov(1,1) = ravecov.dyy();
  cov(1,2) = ravecov.dyz();
  cov(1,3) = ravecov.dypx();
  cov(1,4) = ravecov.dypy();
  cov(1,5) = ravecov.dypz();

  cov(2,0) = ravecov.dxz();
  cov(2,1) = ravecov.dyz();
  cov(2,2) = ravecov.dzz();
  cov(2,3) = ravecov.dzpx();
  cov(2,4) = ravecov.dzpy();
  cov(2,5) = ravecov.dzpz();

  cov(3,0) = ravecov.dxpx();
  cov(3,1) = ravecov.dypx();
  cov(3,2) = ravecov.dzpx();
  cov(3,3) = ravecov.dpxpx();
  cov(3,4) = ravecov.dpxpy();
  cov(3,5) = ravecov.dpxpz();

  cov(4,0) = ravecov.dxpy();
  cov(4,1) = ravecov.dypy();
  cov(4,2) = ravecov.dzpy();
  cov(4,3) = ravecov.dpxpy();
  cov(4,4) = ravecov.dpypy();
  cov(4,5) = ravecov.dpypz();

  cov(5,0) = ravecov.dxpz();
  cov(5,1) = ravecov.dypz();
  cov(5,2) = ravecov.dzpz();
  cov(5,3) = ravecov.dpxpz();
  cov(5,4) = ravecov.dpypz();
  cov(5,5) = ravecov.dpzpz();

  return cov;
}


rave::Point3D
TVector3ToPoint3D(const TVector3 & vec){
  return rave::Point3D(vec.X(), vec.Y(), vec.Z());
}


rave::Covariance3D
TMatrixDSymToCovariance3D(const TMatrixDSym & matrix){
  if (matrix.GetNrows()!=3) {
    Exception exc("TMatrixDSymToCovariance3D ==> TMatrixDSym is not 3x3!",__LINE__,__FILE__);
    throw exc;
  }

  return rave::Covariance3D(matrix(0,0), matrix(0,1), matrix(0,2),
                            matrix(1,1), matrix(1,2), matrix(2,2));

}


} // end of namespace genfit

