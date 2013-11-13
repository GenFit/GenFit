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

#ifndef genfit_mySpacepointDetectorHit_h
#define genfit_mySpacepointDetectorHit_h

#include <TObject.h>


namespace genfit {

/** @brief Example class for a spacepoint detector hit.
 *
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */
class mySpacepointDetectorHit : public TObject {

 public:
  mySpacepointDetectorHit() {;}

  mySpacepointDetectorHit(const TVector3& pos, const TMatrixDSym cov)
  : pos_(pos), cov_(cov) {;}

  const TVector3 getPos() const {return pos_;}
  const TMatrixDSym getCov() const {return cov_;}

 private:

  TVector3 pos_;
  TMatrixDSym cov_;

  ClassDef(mySpacepointDetectorHit,1)
};
/** @} */

} /* End of namespace genfit */

#endif // genfit_mySpacepointDetectorHit_h
