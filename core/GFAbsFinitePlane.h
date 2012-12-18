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
/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */


/** @addtogroup genfit
 * @{
 */

#ifndef GFABSFINITEPLANE_H
#define GFABSFINITEPLANE_H

#include "TObject.h"

/** @brief Abstract base class for implementing arbitrarily shaped finite detector planes
 *
 * Finite detector planes could have any shape, e.g. rectangular, trapezoidal, ...
 * Define them by inheriting from this class. Example, see GFRectFinitePlane
 * This is an optional feature for the GFDetPlane class for avoiding fake intersection points in
 * the fitting of curlers.
 */

class GFAbsFinitePlane : public TObject {
 public:
  //! Returns whether a u,v point is in the active plane or not. Needs to be implemented
  //! in child class.
  virtual bool inActive(const double& u,const double& v)const = 0;

  virtual void Print(const Option_t* = "") const = 0;

  //! Deep copy ctor for polymorphic class.
  virtual GFAbsFinitePlane* clone() const = 0;

  virtual ~GFAbsFinitePlane();

  ClassDef(GFAbsFinitePlane,1)
};

#endif

/** @} */
