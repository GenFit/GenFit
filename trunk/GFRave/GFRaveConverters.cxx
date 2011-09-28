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

#include <iostream>



std::vector < rave::Track >
gfrave::GFTracksToTracks(const std::vector < GFTrack* >  & GFTracks,
                         std::map<unsigned int, GFTrack*>* IdGFTrackMap,
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
  }

  return ravetracks;
}


rave::Track
gfrave::GFTrackToTrack(const GFTrack* orig, int id, std::string tag){
  return gfrave::RepToTrack(orig->getCardinalRep(), id, orig, tag);
}


rave::Track
gfrave::RepToTrack(GFAbsTrackRep* rep, const rave::Track & orig) {
  return gfrave::RepToTrack(rep, orig.id(), orig.originalObject(), orig.tag());
}


rave::Track
gfrave::RepToTrack(GFAbsTrackRep* rep, int id, const void * originaltrack, std::string tag){

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
                  0, tag); // todo: originaltrack pointer not passed!

  return ret;
}


GFDetPlane
gfrave::PlaneToGFDetPlane(const ravesurf::Plane & rplane) {
  return GFDetPlane(gfrave::Point3DToTVector3(rplane.position()),
                    gfrave::Vector3DToTVector3(rplane.normalVector()) );
}


TVector3
gfrave::Point3DToTVector3(const rave::Point3D & v) {
  return TVector3(v.x(), v.y(), v.z());
}

TVector3
gfrave::Vector3DToTVector3(const rave::Vector3D & v) {
  return TVector3(v.x(), v.y(), v.z());
}



