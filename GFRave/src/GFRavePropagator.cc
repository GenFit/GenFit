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


#include "GFRavePropagator.h"
#include "GFRaveConverters.h"
#include "Exception.h"

#include <iostream>


namespace genfit {

GFRavePropagator::GFRavePropagator() :
    IdGFTrackStateMap_(NULL)
{
  //std::cout << "GFRavePropagator::GFRavePropagator() \n";
}
  
GFRavePropagator*
GFRavePropagator::copy() const
{
  //std::cout << "GFRavePropagator::copy() \n";
  return new GFRavePropagator(*this);
}


GFRavePropagator::~GFRavePropagator()
{
  //std::cout << "GFRavePropagator::~GFRavePropagator() \n";
}
    
    
std::pair < rave::Track, double >
GFRavePropagator::to ( const rave::Track & orig,
                       const ravesurf::Cylinder & rcyl ) const
{
  // todo to be implemented!!
  Exception exc("GFRavePropagator::to (cylinder) ==> not yet implemented!",__LINE__,__FILE__);
  throw exc;
}


std::pair < rave::Track, double >
GFRavePropagator::to ( const rave::Track & orig,
                       const ravesurf::Plane & rplane ) const
{
  // will throw Exception if extrapolation does not work
  double path = IdGFTrackStateMap_->at(orig.id()).state_->extrapolateToPlane(PlaneToGFDetPlane(rplane));

  std::pair < rave::Track, double > ret(GFTrackToTrack(IdGFTrackStateMap_->at(orig.id()), orig.id(), orig.tag()), path);
  return ret;
}


rave::Track
GFRavePropagator::closestTo ( const rave::Track & orig,
                              const rave::Point3D & pt, bool transverse ) const
{

  if (transverse){
    Exception exc("GFRavePropagator::closestTo ==> transverse is true, not implemented!",__LINE__,__FILE__);
    throw exc;
  }

  TVector3 point(Point3DToTVector3(pt));
  IdGFTrackStateMap_->at(orig.id()).state_->extrapolateToPoint(point);

  return GFTrackToTrack(IdGFTrackStateMap_->at(orig.id()), orig.id(), orig.tag());
}


void
GFRavePropagator::setIdGFTrackStateMap(std::map < int, genfit::trackAndState > * map){
  //std::cout << "GFRavePropagator::setIdGFMeasuredStateOnPlaneMap() \n";

  IdGFTrackStateMap_ = map;

  if (IdGFTrackStateMap_==NULL) {
    Exception exc("GFRavePropagator::setIdGFMeasuredStateOnPlaneMap ==> map is NULL!",__LINE__,__FILE__);
    throw exc;
  }
  //std::cout<<"IdGFTrackStateMap_: " << (int)IdGFTrackStateMap_ << std::endl;
}

} /* End of namespace genfit */
