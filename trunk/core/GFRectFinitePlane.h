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


#ifndef GFRECTFINITEPLANE_H
#define GFRECTFINITEPLANE_H

#include "GFAbsFinitePlane.h"

/** @brief Concrete implementation of finitie detector plane for rectangles.
 */


class GFRectFinitePlane : public GFAbsFinitePlane {
public:
  //override inActive & Print methods
  bool inActive(const double& u,const double& v) const;
  void Print(const Option_t* = "") const;

  //! give dimensions of finite rectangle: u1,u2,v1,v2
  GFRectFinitePlane(const double&,const double&,const double&,const double&);
  GFRectFinitePlane();

  virtual ~GFRectFinitePlane();

  GFAbsFinitePlane* clone() const {
      return new GFRectFinitePlane(*this);
  }
 private:
  double fUmin,fUmax,fVmin,fVmax;
 public:
  ClassDef(GFRectFinitePlane,1)
};

#endif

/** @} */
