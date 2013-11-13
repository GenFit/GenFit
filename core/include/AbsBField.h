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

#ifndef genfit_AbsBField_h
#define genfit_AbsBField_h

#include <TVector3.h>


namespace genfit {

/** @brief Abstract Interface to magnetic fields in GENFIT
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 */
class AbsBField {
public:

  AbsBField(){;}
  virtual ~AbsBField(){;}

  /**
   * @brief Get the magneticField [kGauss] at position.
   *
   * Override this in your concrete implementation.
   * Provided for compatibility with old genfit.  Use the other interface to avoid
   * unnecessary TVector3 instantiations.
   */
  virtual TVector3 get(const TVector3& position) const = 0;

  /**
   * @brief Get the magneticField [kGauss] at position.
   *
   * Override this in your concrete implementation.
   */
  virtual void get(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) const { const TVector3& B(this->get(TVector3(posX, posY, posZ))); Bx = B.X(); By = B.Y(); Bz = B.Z(); }
 
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsBField_h
