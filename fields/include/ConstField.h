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
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_ConstField_h
#define genfit_ConstField_h

#include "AbsBField.h"


namespace genfit {

/** @brief Constant Magnetic field
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */
class ConstField : public AbsBField {
 public:
  //! define the constant field in this ctor
  ConstField(double b1, double b2, double b3)
    : field_(b1, b2, b3)
  { ; }

  ConstField(const TVector3& field)
    : field_(field)
  { ; }

  //! return value at position
  TVector3 get(const TVector3& pos) const;
  void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const;

 private:
  TVector3 field_;
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_ConstField_h
