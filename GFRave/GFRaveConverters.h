/* Copyright 2008-2009, Technische Universitaet Muenchen,
   Authors: Johannes Rauch

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

/*
 */

/** @addtogroup GFRave
 * @{
 */


#ifndef GFRAVECONVERTERS_H
#define GFRAVECONVERTERS_H

// overwrite visual c stuff
#define RaveDllExport

#include <rave/Track.h>
#include <rave/Plane.h>

#include <GFAbsTrackRep.h>
#include <GFDetPlane.h>
#include <GFTrack.h>

#include <TVector3.h>

#include <iostream>


namespace gfrave{

  rave::Track GFTrackToTrack(GFTrack* orig);
  rave::Track RepToTrack(GFAbsTrackRep* rep, const rave::Track & orig);
  rave::Track RepToTrack(GFAbsTrackRep* rep, int id, void * originaltrack = 0, std::string tag="");

  GFDetPlane PlaneToGFDetPlane(const ravesurf::Plane & rplane);

  TVector3 Point3DToTVector3(const rave::Point3D &);
  TVector3 Vector3DToTVector3(const rave::Vector3D &);

}


#endif

/** @} */
