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

#include <GFTrack.h>
#include <GFAbsTrackRep.h>

#include <rave/Plane.h>
#include <rave/Point3D.h>

#include <iostream>


rave::Track RepToTrack(GFAbsTrackRep* rep, const rave::Track & orig) {
  return gfrave::RepToTrack(rep, orig->id(), orig->originalObject(), orig->tag());
}


rave::Track RepToTrack(GFAbsTrackRep* rep, int id, void * originaltrack, std::string tag){

  GFDetPlane* refPlane(rep->getReferencePlane());
  TVector3 pos, mom;
  TMatrixT<double> cov;

  rep->getPosMomCov(refPlane, pos, mom, cov);

  // state
  rave::Vector6D state(pos.X(), pos.Y(), pos.Z(),
                       mom.X(), mom.Y(), mom.Z());

  // covariance
  rave::Covariance6D cov(cov[0][0], cov[1][0], cov[2][0],
                         cov[1][1], cov[2][1], cov[2][2],
                         cov[3][0], cov[4][0], cov[5][0],
                         cov[3][1], cov[4][1], cov[5][1],
                         cov[3][2], cov[4][2], cov[5][2],
                         cov[3][3], cov[4][3], cov[5][3],
                         cov[4][4], cov[5][4], cov[5][5]);

  rave::Track ret(id, state, cov,
                  rep->getCharge(), rep->getChiSqu(), rep->getNDF(),
                  originaltrack, tag);

  return ret;
}


GFDetPlane PlaneToGFDetPlane(const ravesurf::Plane & rplane) {
  return GFDetPlane(Point3DToTVector3(rplane.position()),
                    Vector3DToTVector3(rplane.normalVector()) );
}


TVector3 Point3DToTVector3(const rave::Point3D & v) {
  return TVector3(v.x(), v.y(), v.z());
}

TVector3 Vector3DToTVector3(const rave::Vector3D & v) {
  return TVector3(v.x(), v.y(), v.z());
}
