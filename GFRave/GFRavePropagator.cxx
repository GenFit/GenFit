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
#include "GFRaveConverters.h"
#include <GFException.h>

#include <iostream>

using namespace std;


gfrave::GFRavePropagator::GFRavePropagator() :
  IdGFTrackMap(NULL)
{}
  
gfrave::GFRavePropagator*
gfrave::GFRavePropagator::copy() const
{
  return new gfrave::GFRavePropagator ( *this );
}


gfrave::GFRavePropagator::~GFRavePropagator()
{

}
    
    
std::pair < rave::Track, double >
gfrave::GFRavePropagator::to ( const rave::Track & orig,
                               const ravesurf::Cylinder & rcyl ) const
{
  // todo to be implemented!!
}


std::pair < rave::Track, double >
gfrave::GFRavePropagator::to ( const rave::Track & orig,
                               const ravesurf::Plane & rplane ) const
{
  GFAbsTrackRep* rep = checkTrack(orig);
  double path = rep->extrapolate(gfrave::PlaneToGFDetPlane(rplane));

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::to ==> Extrapolation failed!",__LINE__,__FILE__);
    throw exc;
  }
  
  std::pair < rave::Track, double > ret(gfrave::RepToTrack(rep, orig), path);
  return ret;
}


rave::Track
gfrave::GFRavePropagator::closestTo ( const rave::Track & orig,
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




GFAbsTrackRep*
gfrave::GFRavePropagator::checkTrack(const rave::Track & track) const{
  if (IdGFTrackMap==NULL) {
    GFException exc("GFRavePropagator::checkTrack ==> IdGFTrackMap is not defined, cannot access GFTracks!",__LINE__,__FILE__);
    throw exc;
  }

  if (!(track.isValid())) {
    GFException exc("GFRavePropagator::checkTrack ==> Track is not valid!",__LINE__,__FILE__);
    throw exc;
  }

  if (IdGFTrackMap->count(track.id()) == 0) {
    GFException exc("GFRavePropagator::checkTrack ==> no entry in IdGFTrackMap corresponding to track id, cannot access corresponding GFTrack!",__LINE__,__FILE__);
    throw exc;
  }

  GFAbsTrackRep* rep = IdGFTrackMap->at(track.id())->getCardinalRep();

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::checkTrack ==> Status flag != 0, cannot extrapolate!",__LINE__,__FILE__);
    throw exc;
  }

  return rep;
}

