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
/** @addtogroup genfit
 * @{
 */

#ifndef GFABSPROLATESPACEPOINT_H
#define GFABSPROLATESPACEPOINT_H


#include "GFAbsSpacepointHit.h"

/** @brief Abstract hit class implementing a space point hit geometry with a very prolate
 * form of the covariance matrix.
 *
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * RecoHits for detectors measuring 3D space points with errors in one direction
 * much larger than the errors perpendicular should inherit from GFAbsProlateSpacepointHit.
 *
 * For these hits, a virtual detector plane lying in the POCA and
 * perpendicular to the track yields wrong results. Instead, the plane should contain the
 * direction of the largest error.
 *
 * The largest error direction can be set. Standard is in z.
 *
 */

class GFAbsProlateSpacepointHit : public GFAbsSpacepointHit {
 public:

  // Constructors/Destructors ---------
  GFAbsProlateSpacepointHit(unsigned int dim = 3);

  virtual ~GFAbsProlateSpacepointHit(){;}
  
  // Operations ----------------------
   /** @brief Get detector plane perpendicular to track.
    *
    * The detector plane is contructed from the position of the hit and
    * the track representation. For this the track is extrapolated to the
    * point of closest approach to a line parallel to the largest error direction.
    */
  virtual const GFDetPlane& getDetPlane(GFAbsTrackRep* rep);

  const TVector3& getLargestErrorDirection(){return fLargestErrorDirection;}
  void setLargestErrorDirection(TVector3& dir){fLargestErrorDirection = dir.Unit();}

 protected:

  // Private Data Members ------------
  TVector3 fLargestErrorDirection; // direction of largest error

  // Private Methods -----------------

 public:
  ClassDef(GFAbsProlateSpacepointHit,1);

};

#endif

/** @} */
