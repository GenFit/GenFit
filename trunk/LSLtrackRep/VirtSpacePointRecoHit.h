/* Copyright 2008-2009, Technische Universitaet Muenchen,
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
 *  @{ */

/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */
#ifndef VIRTSPACEPOINTRECOHIT_H
#define VIRTSPACEPOINTRECOHIT_H

#include "GFRecoHitIfc.h"
#include "GFSpacepointHitPolicy.h"

/** @brief Define a hit policy template specialization */
typedef GFRecoHitIfc<GFSpacepointHitPolicy> SpacepointRecoHit;

/** @brief A basic space point reco hit. Uses a TVector3 for initialization
 *
 * This class can be used to calculate the distance (residual) of a track
 * to a point in space.
 *
 * It is a very basic example of a RecoHit. VirtSpacePointRecoHit is 
 * initialized by an ordinary 3-vector and has unity error matrix. 
 * It is also a good example which can be used to develop more 
 * realistic RecoHits.
 */
class VirtSpacePointRecoHit : public SpacepointRecoHit {
public:
 
  // Constructors/Destructors ---------
  VirtSpacePointRecoHit();
  VirtSpacePointRecoHit(double x, double y, double z);
  /** @brief Initializing constructor
   *
   * With a RecoHitProducer that works on a TClonseArray containing 
   * TVector3 objects this constructor allows the automatic conversion
   * from a cluster class (here the TVector3) to a RecoHit 
   * (the VirtSpacePointRecoHit). 
   */
  VirtSpacePointRecoHit(const TVector3& pos);

  virtual ~VirtSpacePointRecoHit();

  virtual GFAbsRecoHit* clone();

  // Operations ----------------------
  /** @brief currently works only with LSLTrackRep
   *
   * Modifications needed here to use VirtSpacePointRecoHit with different
   * track representations.
   *
   * @sa GFAbsRecoHit::getHMatrix
   */
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);


private:

  // Private Data Members ------------
  static const int NparHitRep = 3;

  // Private Methods -----------------

public:
  ClassDef(VirtSpacePointRecoHit,1)


};

#endif

/** @} */
