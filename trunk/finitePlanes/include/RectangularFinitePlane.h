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

#ifndef genfit_RectangularFinitePlane_h
#define genfit_RectangularFinitePlane_h

#include "AbsFinitePlane.h"


namespace genfit {

/**
 * @brief Rectangular finite plane.
 */
class RectangularFinitePlane : public AbsFinitePlane {

 public:

  //! give dimensions of finite rectangle: u1,u2, v1,v2
  RectangularFinitePlane(const double&, const double&, const double&, const double&);
  RectangularFinitePlane();
  virtual ~RectangularFinitePlane();

  //override inActive & Print methods
  bool isInActive(double u, double v) const;
  void Print(const Option_t* = "") const;

  RectangularFinitePlane* clone() const {
    return new RectangularFinitePlane(*this);
  }

 private:

  double uMin_, uMax_, vMin_, vMax_;

 public:

  ClassDef(RectangularFinitePlane,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_RectangularFinitePlane_h
