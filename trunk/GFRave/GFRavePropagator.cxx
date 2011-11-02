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
#include "GFException.h"

#include <iostream>

using namespace std;


GFRavePropagator::GFRavePropagator() :
    IdGFTrackRepMap(NULL)
{}
  
GFRavePropagator*
GFRavePropagator::copy() const
{
  return new GFRavePropagator(*this);
}


GFRavePropagator::~GFRavePropagator()
{

}
    
    
std::pair < rave::Track, double >
GFRavePropagator::to ( const rave::Track & orig,
                       const ravesurf::Cylinder & rcyl ) const
{
  // todo to be implemented!!
  GFException exc("GFRavePropagator::to (cylinder) ==> not yet implemented!",__LINE__,__FILE__);
  throw exc;
}


std::pair < rave::Track, double >
GFRavePropagator::to ( const rave::Track & orig,
                       const ravesurf::Plane & rplane ) const
{
  GFAbsTrackRep* rep = getTrackRep(orig);
  double path = rep->extrapolate(GFRave::PlaneToGFDetPlane(rplane));

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::to ==> Extrapolation failed!",__LINE__,__FILE__);
    throw exc;
  }
  
  std::pair < rave::Track, double > ret(GFRave::RepToTrack(rep, orig), path);
  return ret;
}


rave::Track
GFRavePropagator::closestTo ( const rave::Track & orig,
                              const rave::Point3D & pt, bool transverse ) const
{

  if (transverse){
    GFException exc("GFRavePropagator::closestTo ==> transverse is true, not implemented!",__LINE__,__FILE__);
    throw exc;
  }

  GFAbsTrackRep* rep = getTrackRep(orig);

  TVector3 poca, normvec;
  TVector3 point(GFRave::Point3DToTVector3(pt));
  rep->extrapolateToPoint(point, poca, normvec);

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::to ==> Extrapolation failed!",__LINE__,__FILE__);
    throw exc;
  }

  return GFRave::RepToTrack(rep, orig);
}


GFAbsTrackRep*
GFRavePropagator::getTrackRep(const rave::Track & track) const{

  if (IdGFTrackRepMap==NULL) {
    GFException exc("GFRavePropagator::getTrackRep ==> IdGFTrackRepMap is NULL, cannot access GFTracks!",__LINE__,__FILE__);
    throw exc;
  }

  if (!(track.isValid())) {
    GFException exc("GFRavePropagator::getTrackRep ==> Track is not valid!",__LINE__,__FILE__);
    throw exc;
  }

  std::cerr<<"GFRavePropagator::getTrackRep track id: "<<track.id()<<std::endl;
  std::cerr<<"  pos: "; GFRave::Point3DToTVector3(track.state().position()).Print();
  std::cerr<<"  mom: "; GFRave::Vector3DToTVector3(track.state().momentum()).Print();

  if (IdGFTrackRepMap->count(track.id()) == 0) {
    GFException exc("GFRavePropagator::getTrackRep ==> no entry in IdGFTrackRepMap corresponding to track id, cannot access corresponding track representation!",__LINE__,__FILE__);
    throw exc;
  }

  GFAbsTrackRep* rep = IdGFTrackRepMap->at(track.id());

  if (rep->getStatusFlag() != 0){
    GFException exc("GFRavePropagator::getTrackRep ==> Status flag != 0, cannot extrapolate!",__LINE__,__FILE__);
    throw exc;
  }

  GFRave::setTrackRepData(track, rep); // set trackrep state and cov

  rep->setStatusFlag(0);

  return rep;
}


void
GFRavePropagator::setIdGFTrackRepMap(std::map < int, GFAbsTrackRep* > * map){
  if (map==NULL) {
    GFException exc("GFRavePropagator::setIdGFTrackRepMap ==> map is NULL!",__LINE__,__FILE__);
    throw exc;
  }
  IdGFTrackRepMap = map;
  if (IdGFTrackRepMap==NULL) {
    GFException exc("GFRavePropagator::setIdGFTrackRepMap ==> IdGFTrackMap is NULL!",__LINE__,__FILE__);
    throw exc;
  }
  std::cout<<"IdGFTrackRepMap: " << (int)IdGFTrackRepMap << std::endl;
}

