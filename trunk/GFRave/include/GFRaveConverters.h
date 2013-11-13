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

/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup GFRave
 * @{
 */


#ifndef GFRAVECONVERTERS_H
#define GFRAVECONVERTERS_H

#include "AbsTrackRep.h"
#include "DetPlane.h"
#include "Track.h"
#include "GFRaveVertex.h"
#include "GFRaveVertexFactory.h"

#include <rave/Track.h>
#include <rave/Plane.h>
#include <rave/Vertex.h>

#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>

#include <iostream>


/**
 * @brief Collection of converter functions
 **/
namespace genfit {

  // GENFIT to RAVE

  /** Convert a vector of genfit::Tracks to rave::Tracks
   * also builds a map of unique ids to genfit::Tracks; These ids are stored in the rave::tracks. They are counted from startID.
   * The map has to be passed to the GFRavePropagator, so that it can access the trackreps of the genfit::Tracks corresponding to the
   * rave::tracks.
   */
  std::vector < rave::Track > GFTracksToTracks(const std::vector < genfit::Track* > & GFTracks,
                                               std::map<int, genfit::trackAndState>& IdGFTrackStateMap,
                                               int startID = 0);

  rave::Track GFTrackToTrack(trackAndState, int id = -1, std::string tag="");
  //rave::Track MeasuredStateOnPlaneToTrack(const MeasuredStateOnPlane* state, const rave::Track& orig);
  //rave::Track MeasuredStateOnPlaneToTrack(const MeasuredStateOnPlane* state, int id = -1, Track* originaltrack = NULL, std::string tag="");

  // RAVE to GENFIT
  /** set state and cov of a MeasuredStateOnPlane according to rave track
   */
  void setData(const rave::Track & orig, MeasuredStateOnPlane* state);

  GFRaveVertex* RaveToGFVertex(const rave::Vertex &, const std::map<int, genfit::trackAndState>& IdGFTrackStateMap);
  void RaveToGFVertices(std::vector<GFRaveVertex*> *, const std::vector<rave::Vertex> &, const std::map<int, genfit::trackAndState>& IdGFTrackStateMap);

  SharedPlanePtr PlaneToGFDetPlane(const ravesurf::Plane& rplane);

  // RAVE to ROOT
  TVector3 Point3DToTVector3(const rave::Point3D&);
  TVector3 Vector3DToTVector3(const rave::Vector3D&);

  TMatrixDSym Covariance3DToTMatrixDSym(const rave::Covariance3D&);
  TVectorD Vector6DToTVectorD(const rave::Vector6D&);
  TMatrixDSym Covariance6DToTMatrixDSym(const rave::Covariance6D&);

  // ROOT to RAVE
  rave::Point3D TVector3ToPoint3D(const TVector3 &);
  rave::Covariance3D TMatrixDSymToCovariance3D(const TMatrixDSym&);


} /* End of namespace genfit */
/** @} */

#endif // GFRAVECONVERTERS_H

