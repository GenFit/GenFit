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


#include "GFRaveConverters.h"

#include <GFTrack.h>
#include <GFAbsTrackRep.h>
#include <GFException.h>

#include <rave/Plane.h>

#include "GFRaveTrackParameters.h"

#include <iostream>



std::vector < rave::Track >
GFRave::GFTracksToTracks(const std::vector < GFTrack* >  & GFTracks,
                         std::map<int, GFTrack*>* IdGFTrackMap,
                         std::map<int, GFAbsTrackRep*> * IdGFTrackRepMap,
                         int startID){

  unsigned int ntracks(GFTracks.size());

  std::vector < rave::Track > ravetracks;
  ravetracks.reserve(ntracks);

  for (unsigned int i=0; i<ntracks; ++i){
    ravetracks.push_back(GFTrackToTrack(GFTracks[i], startID++) );

    if (IdGFTrackMap != NULL){
      if (IdGFTrackMap->count(startID)>0){
        GFException exc("GFTracksToTracks ==> IdGFTrackMap has already an entry for this id",__LINE__,__FILE__);
        throw exc;
      }
      (*IdGFTrackMap)[startID] = GFTracks[i];
    }

    if (IdGFTrackRepMap != NULL){
      if (IdGFTrackRepMap->count(startID)>0){
        GFException exc("GFTracksToTracks ==> IdGFTrackRepMap has already an entry for this id",__LINE__,__FILE__);
        throw exc;
      }
      (*IdGFTrackRepMap)[startID] = GFTracks[i]->getCardinalRep()->clone(); // here clones are made so that the state of the original GFTracks and their TrackReps will not be altered by the vertexing process
    }

  }

  return ravetracks;
}


rave::Track
GFRave::GFTrackToTrack(GFTrack* orig, int id, std::string tag){
  return GFRave::RepToTrack(orig->getCardinalRep(), id, orig, tag);
}


rave::Track
GFRave::RepToTrack(GFAbsTrackRep* rep, const rave::Track & orig) {
  return GFRave::RepToTrack(rep, orig.id(), orig.originalObject(), orig.tag());
}


rave::Track
GFRave::RepToTrack(GFAbsTrackRep* rep, int id, void * originaltrack, std::string tag){

  GFDetPlane refPlane(rep->getReferencePlane());
  TVector3 pos, mom;
  TMatrixT<double> cov;

  rep->getPosMomCov(refPlane, pos, mom, cov);

  // state
  rave::Vector6D ravestate(pos.X(), pos.Y(), pos.Z(),
                           mom.X(), mom.Y(), mom.Z());

  // covariance
  rave::Covariance6D ravecov(cov[0][0], cov[1][0], cov[2][0],
                             cov[1][1], cov[2][1], cov[2][2],
                             cov[3][0], cov[4][0], cov[5][0],
                             cov[3][1], cov[4][1], cov[5][1],
                             cov[3][2], cov[4][2], cov[5][2],
                             cov[3][3], cov[4][3], cov[5][3],
                             cov[4][4], cov[5][4], cov[5][5]);

  rave::Track ret(id, ravestate, ravecov,
                  rep->getCharge(), rep->getChiSqu(), rep->getNDF(),
                  originaltrack, tag);

  return ret;
}


void
GFRave::setTrackRepData(const rave::Track & orig, GFAbsTrackRep* rep){

  rep->setPosMomCov(TVector3(orig.state().x(), orig.state().y(), orig.state().z()),
                    TVector3(orig.state().px(), orig.state().py(), orig.state().pz()),
                    Covariance6DToTMatrixT(orig.error()));

}


GFRaveVertex*
GFRave::RaveToGFVertex(const rave::Vertex & ravevertex, const std::map<int, GFTrack*> * IdGFTrackMap){

  if (!(ravevertex.isValid())) {
    GFException exc("RaveToGFVertex ==> rave Vertex is not valid!",__LINE__,__FILE__);
    throw exc;
  }

  std::vector < std::pair < double, GFTrack* > > originalTracks;
  std::vector < std::pair < double, GFRaveTrackParameters > > smoothedTracks;

  if (IdGFTrackMap!=NULL){
    std::vector < std::pair < float, rave::Track > > raveweightedtracks = ravevertex.weightedTracks();
    std::vector < std::pair < float, rave::Track > > ravesmoothedtracks = ravevertex.weightedRefittedTracks();

    for (unsigned int i=0; i<raveweightedtracks.size(); ++i){
      GFTrack* origtrack = IdGFTrackMap->at(raveweightedtracks[i].second.id());
      originalTracks.push_back(std::pair< double, GFTrack* > (raveweightedtracks[i].first, origtrack) );

    }

    for (unsigned int i=0; i<ravesmoothedtracks.size(); ++i){
      GFTrack* origtrack = IdGFTrackMap->at(raveweightedtracks[i].second.id());

      GFRaveTrackParameters trackparams(GFRave::Vector6DToTMatrixT(raveweightedtracks[i].second.state()),
                                        GFRave::Covariance6DToTMatrixT(raveweightedtracks[i].second.error()),
                                        origtrack->getCardinalRep()->getCharge(),
                                        origtrack->getCardinalRep()->getPDG());

      smoothedTracks.push_back(std::pair < double, GFRaveTrackParameters > (ravesmoothedtracks[i].first, trackparams) );

    }

  }

  return new GFRaveVertex(GFRave::Point3DToTVector3(ravevertex.position()),
                          GFRave::Covariance3DToTMatrixT(ravevertex.error()),
                          originalTracks,
                          smoothedTracks,
                          ravevertex.ndf(), ravevertex.chiSquared(), ravevertex.id());
}


std::vector<GFRaveVertex*> *
GFRave::RaveToGFVertices(const std::vector<rave::Vertex> & raveVertices, const std::map<int, GFTrack*> * IdGFTrackMap){

  unsigned int nVert(raveVertices.size());

  std::vector<GFRaveVertex*> * GFVertices = new std::vector<GFRaveVertex*>;
  GFVertices->reserve(nVert);

  for (unsigned int i=0; i<nVert; ++i){
    GFVertices->push_back(RaveToGFVertex(raveVertices[i], IdGFTrackMap));
  }

  return GFVertices;
}


GFDetPlane
GFRave::PlaneToGFDetPlane(const ravesurf::Plane & rplane) {
  return GFDetPlane(GFRave::Point3DToTVector3(rplane.position()),
                    GFRave::Vector3DToTVector3(rplane.normalVector()) );
}


TVector3
GFRave::Point3DToTVector3(const rave::Point3D & v) {
  return TVector3(v.x(), v.y(), v.z());
}

TVector3
GFRave::Vector3DToTVector3(const rave::Vector3D & v) {
  return TVector3(v.x(), v.y(), v.z());
}


TMatrixT<double>
GFRave::Covariance3DToTMatrixT(const rave::Covariance3D & ravecov){
  TMatrixT<double> cov(3,3);

  cov[0][0] = ravecov.dxx();
  cov[1][1] = ravecov.dyy();
  cov[2][2] = ravecov.dzz();

  cov[0][1] = ravecov.dxy();
  cov[1][0] = ravecov.dxy();

  cov[0][2] = ravecov.dxz();
  cov[2][0] = ravecov.dxz();

  cov[1][2] = ravecov.dyz();
  cov[2][1] = ravecov.dyz();

  assert (cov.IsSymmetric()==true); // todo: for QA reasons

  return cov;
}


TMatrixT<double>
GFRave::Vector6DToTMatrixT(const rave::Vector6D & ravevec){
  TMatrixT<double> vec(1,6);

  vec[0][0] = ravevec.x();
  vec[0][1] = ravevec.y();
  vec[0][2] = ravevec.z();

  vec[0][3] = ravevec.px();
  vec[0][4] = ravevec.py();
  vec[0][5] = ravevec.pz();

  return vec;
}


TMatrixT<double>
GFRave::Covariance6DToTMatrixT(const rave::Covariance6D & ravecov){
  TMatrixT<double> cov(6,6);

  cov[0][0] = ravecov.dxx();
  cov[1][1] = ravecov.dyy();
  cov[2][2] = ravecov.dzz();

  cov[0][1] = ravecov.dxy();
  cov[1][0] = ravecov.dxy();

  cov[0][2] = ravecov.dxz();
  cov[2][0] = ravecov.dxz();

  cov[1][2] = ravecov.dyz();
  cov[2][1] = ravecov.dyz();


  cov[3][3] = ravecov.dpxpx();
  cov[4][4] = ravecov.dpypy();
  cov[5][5] = ravecov.dpzpz();

  cov[3][4] = ravecov.dpxpy();
  cov[4][3] = ravecov.dpxpy();

  cov[3][5] = ravecov.dpxpz();
  cov[5][3] = ravecov.dpxpz();

  cov[4][5] = ravecov.dpypz();
  cov[5][4] = ravecov.dpypz();



  cov[3][0] = ravecov.dxpx();
  cov[0][3] = ravecov.dxpx();
  cov[4][1] = ravecov.dypy();
  cov[1][4] = ravecov.dypy();
  cov[5][2] = ravecov.dzpz();
  cov[2][5] = ravecov.dzpz();

  cov[4][0] = ravecov.dxpy();
  cov[0][4] = ravecov.dxpy();

  cov[5][0] = ravecov.dxpz();
  cov[0][5] = ravecov.dxpz();

  cov[3][1] = ravecov.dypx();
  cov[1][3] = ravecov.dypx();

  cov[5][1] = ravecov.dypz();
  cov[1][5] = ravecov.dypz();

  cov[4][2] = ravecov.dzpy();
  cov[2][4] = ravecov.dzpy();

  cov[5][2] = ravecov.dzpz();
  cov[2][5] = ravecov.dzpz();

  assert (cov.IsSymmetric()==true); // todo: for QA reasons

  return cov;
}


