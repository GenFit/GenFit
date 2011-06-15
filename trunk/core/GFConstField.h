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

#ifndef GFCONSTFIELD_H
#define GFCONSTFIELD_H

#include"GFAbsBField.h"

/** @brief Constant Magnetic field
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */

class GFConstField : public GFAbsBField{
 public:
  //! define the constant field in this ctor
  GFConstField(double b1,double b2, double b3){
    fF1 = b1;    fF2 = b2;    fF3 = b3;
  }

  //! return value at position
  TVector3 get(const TVector3& pos) const;

 private:
  double fF1,fF2,fF3;
};

/** @} */ 
#endif
