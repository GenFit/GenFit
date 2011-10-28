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

#ifndef GFPSEUDOSPACEPOINTWIREHITPOLICY_H
#define GFPSEUDOSPACEPOINTWIREHITPOLICY_H


#include "TMatrixT.h"
#include "TObject.h"

#include "GFDetPlane.h"

class GFAbsRecoHit;
class GFAbsTrackRep;

/** @brief Policy class implementing a poseudo space point hit geometry.
 *
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * RecoHits for detectors measuring 3D space points with errors in one direction (the "wire direction")
 * much larger than the errors perpendicular should inherit from RecoHitIfc<GFPseudoSpacepointWireHitPolicy>.
 *
 * For these hits, a virtual detector plane lying in the POCA and
 * perpendicular to the track yields wrong results. Instead, the plane should contain the
 * direction of the largest error (i.e. the "wire direction").
 *
 * The "wire direction" can be set. Standard is in z.
 *
 */

class GFPseudoSpacepointWireHitPolicy : public TObject {
public:


  GFPseudoSpacepointWireHitPolicy();
  
  // Operations ----------------------
   /** @brief Get detector plane perpendicular to track.
    *
    * The detector plane is contructed from the position of the hit and
    * the track representation. For this the track is extrapolated to the
    * point of closest approach to a line parallel to z and containing the hit.
    */
  const GFDetPlane& detPlane(GFAbsRecoHit*, GFAbsTrackRep*);

  /** @brief Hit coordinates in detector plane.
   */
  TMatrixT<double> hitCoord(GFAbsRecoHit*,const GFDetPlane&);

  /** @brief Hit covariances in detector plane.
   */
  TMatrixT<double> hitCov(GFAbsRecoHit*,const GFDetPlane&);

  virtual ~GFPseudoSpacepointWireHitPolicy(){;}

  double getMaxDistance(){return fMaxdistance;}
  void setMaxDistance(double d){fMaxdistance=d;}

  TVector3 getWireDirection(){return fWireDirection;}
  void setWireDirection(TVector3 & dir){fWireDirection = dir; fWireDirection.SetMag(1.);}

  const std::string& getName(){return fPolicyName;}

 private:
  static const std::string fPolicyName;

  // Private Data Members ------------
  GFDetPlane fDetPlane;
  double fMaxdistance;
  TVector3 fWireDirection; // direction of largest error

  // Private Methods -----------------

 public:
  ClassDef(GFPseudoSpacepointWireHitPolicy,1);
};

#endif

/** @} */
