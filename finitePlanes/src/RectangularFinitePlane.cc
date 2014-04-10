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

#include <iostream>
#include <cassert>

#include "RectangularFinitePlane.h"

namespace genfit {

RectangularFinitePlane::RectangularFinitePlane(const double& umin,const double& umax,
             const double& vmin,const double& vmax)
  : uMin_(umin),uMax_(umax),vMin_(vmin),vMax_(vmax)
{
  assert(umin<umax);
  assert(vmin<vmax);
}

RectangularFinitePlane::RectangularFinitePlane()
  : uMin_(1.),uMax_(-1.),vMin_(1.),vMax_(-1.)//for this default ctor inActive always false
{}


RectangularFinitePlane::~RectangularFinitePlane(){

}

bool RectangularFinitePlane::isInActive(double u, double v) const{
  return (u>=uMin_ && u<=uMax_ && v>=vMin_ && v<=vMax_);
}

void RectangularFinitePlane::Print(const Option_t*) const{
  std::cout << "Rectangular Finite Plane Umin=" << uMin_ << ", Umax="
      << uMax_ << ", Vmin=" << vMin_ << ", Vmax=" << vMax_ << std::endl;
}

} /* End of namespace genfit */
