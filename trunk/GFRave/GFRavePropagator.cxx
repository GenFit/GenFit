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


#include "GFRavePropagator.h"

#include <iostream>

using namespace std;


rave::GFRavePropagator::GFRavePropagator() :
  IdGFTrackMap(NULL)
{}
  
rave::Propagator * rave::GFRavePropagator::copy() const
{
  return new rave::GFRavePropagator ( * this );
}


rave::GFRavePropagator::~GFRavePropagator()
{

}
    
    
rave::Track rave::GFRavePropagator::to ( const rave::Track & orig,
                          const ravesurf::Cylinder & rcyl ) const
{
}


rave::Track rave::GFRavePropagator::to ( const rave::Track & orig,
                          const ravesurf::Plane & rplane ) const
{
  GFAbsTrackRep* rep = checkTrack(orig);
  rep->extrapolateToPlane(gfrave::PlaneToGFDetPlane(rplane));

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::to ==> Extrapolation failed!",__LINE__,__FILE__);
    throw exc;
  }
  
  return gfrave::RepToTrack(rep, orig);
}


rave::Track rave::GFRavePropagator::closestTo ( const rave::Track & orig,
    const rave::Point3D & pt, bool transverse ) const
{
  GFAbsTrackRep* rep = checkTrack(orig);

  TVector3 poca, normvec;
  rep->extrapolateToPoint(gfrave::Point3DToTVector3(pt), poca, normvec);

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::to ==> Extrapolation failed!",__LINE__,__FILE__);
    throw exc;
  }

  return gfrave::RepToTrack(rep, orig);
}




GFAbsTrackRep* checkTrack(const rave::Track & track) const{
  if (IdGFTrackMap==NULL) {
    GFException exc("GFRavePropagator::checkTrack ==> IdGFTrackMap is not defined, cannot access GFTracks!",__LINE__,__FILE__);
    throw exc;
  }

  if (!(track->isValid())) {
    GFException exc("GFRavePropagator::checkTrack ==> Track is not valid!",__LINE__,__FILE__);
    throw exc;
  }

  if (IdGFTrackMap->count(track->id()) == 0) {
    GFException exc("GFRavePropagator::checkTrack ==> no entry in IdGFTrackMap corresponding to track id, cannot access corresponding GFTrack!",__LINE__,__FILE__);
    throw exc;
  }

  rep = (*IdGFTrackMap)[track->id()];

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::checkTrack ==> Status flag != 0, cannot extrapolate!",__LINE__,__FILE__);
    throw exc;
  }

  return rep;
}

