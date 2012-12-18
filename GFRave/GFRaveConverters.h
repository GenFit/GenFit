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

#include "rave/Track.h"
#include "rave/Plane.h"
#include "rave/Vertex.h"

#include "GFAbsTrackRep.h"
#include "GFDetPlane.h"
#include "GFTrack.h"
#include "GFRaveVertex.h"

#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include <iostream>

/**
 * @brief Collection of converter functions
 **/

namespace GFRave
{


  // GENFIT to RAVE

  /** Convert a vector of GFTracks to rave::Tracks
   * also builds a map of unique ids to GFTracks; These ids are stored in the rave::tracks. They are counted from startID.
   * The map has to be passed to the GFRavePropagator, so that it can access the trackreps of the GFTracks corresponding to the
   * rave::tracks.
   */
  std::vector < rave::Track > GFTracksToTracks(const std::vector < GFTrack* > & GFTracks,
                                               std::map<int, GFTrack*> * IdGFTrackMap = NULL,
                                               std::map<int, GFAbsTrackRep*> * IdGFTrackRepMap = NULL,
                                               int startID = 0);

  std::vector < rave::Track > GFTrackRepsToTracks(const std::vector < GFAbsTrackRep* >  & GFTrackReps,
                                                  std::map<int, GFTrack*> * IdGFTrackMap,
                                                  std::map<int, GFAbsTrackRep*> * IdGFTrackRepMap,
                                                  int startID);

  rave::Track GFTrackToTrack(GFTrack* orig, int id = -1, std::string tag="");
  rave::Track RepToTrack(GFAbsTrackRep* rep, const rave::Track & orig);
  rave::Track RepToTrack(GFAbsTrackRep* rep, int id = -1, void * originaltrack = 0, std::string tag="");

  // RAVE to GENFIT
  /** sets state and cov of a trackrep according to rave track
   */
  void setTrackRepData(const rave::Track & orig, GFAbsTrackRep* rep);

  GFRaveVertex* RaveToGFVertex(const rave::Vertex &, const std::map<int, GFTrack*> * IdGFTrackMap,  const std::map<int, GFAbsTrackRep*> * IdGFTrackRepMap);
  void RaveToGFVertices(std::vector<GFRaveVertex*> *, const std::vector<rave::Vertex> &, const std::map<int, GFTrack*> * IdGFTrackMap,  const std::map<int, GFAbsTrackRep*> * IdGFTrackRepMap);

  GFDetPlane PlaneToGFDetPlane(const ravesurf::Plane & rplane);

  // RAVE to ROOT
  TVector3 Point3DToTVector3(const rave::Point3D &);
  TVector3 Vector3DToTVector3(const rave::Vector3D &);

  TMatrixDSym Covariance3DToTMatrixDSym(const rave::Covariance3D &);
  TVectorD Vector6DToTVectorD(const rave::Vector6D &);
  TMatrixDSym Covariance6DToTMatrixDSym(const rave::Covariance6D &);

  // ROOT to RAVE
  rave::Point3D TVector3ToPoint3D(const TVector3 &);
  rave::Covariance3D TMatrixDSymToCovariance3D(const TMatrixDSym &);

}


#endif

/** @} */
