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

#ifndef GFDAFWIREHIT_H
#define GFDAFWIREHIT_H

#include "GFDafHit.h"
#include "GFAbsWireHit.h"

/** @brief Wrapper class for use with GFDaf.
 *
 * The GFDafWireHit is a hit class which acts as an effective hit.
 * It basically makes two (a left and a right) hits out of a wire hit,
 * so that the DAF can resolve the left/right ambiguity.
 * Initial weights will be set if according to GFAbsWireHit::getLeftRightResolution().
 */

class GFDafWireHit: public GFDafHit {
 public:

	/** @brief Constructor taking a vector of ONE wire hit.
	 */
  GFDafWireHit(GFAbsWireHit* hit);

  ~GFDafWireHit();

	/** @brief Get the measurement m,V
	 */
  void getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V);

  virtual void getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V, unsigned int iHit);

 public:
	ClassDef(GFDafWireHit,3)

};
#endif

/** @} */

