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
/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */


/** @addtogroup genfit
 * @{
 */

#ifndef genfit_AbsFinitePlane_h
#define genfit_AbsFinitePlane_h

#include <Rtypes.h> // for Option_t


namespace genfit {

/**
 * @brief Abstract base class for finite detector planes.
 *
 * This is most important for avoiding fake intersection points in fitting of curlers.
 * This should be implemented for silicon detectors most importantly.
 */
class AbsFinitePlane {

 public:

  AbsFinitePlane() {};
  virtual ~AbsFinitePlane() {};

  //! Returns whether a u,v point is in the active plane or not. Needs to be implemented
  //! in child class.
  virtual bool isInActive(double u, double v) const = 0;

  //! Deep copy ctor for polymorphic class.
  virtual AbsFinitePlane* clone() const = 0;

  virtual void Print(const Option_t* = "") const = 0;


 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsFinitePlane(const AbsFinitePlane&) {;}
  AbsFinitePlane& operator=(const AbsFinitePlane&);
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsFinitePlane_h
