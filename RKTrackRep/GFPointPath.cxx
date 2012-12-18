/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "GFPointPath.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>

GFPointPath::GFPointPath()
  : fX(0),  fY(0),  fZ(0),
    fPath(0)
{
  ;
}

GFPointPath::GFPointPath(const TVector3& pos, double path)
  : fX(pos.X()),  fY(pos.Y()),  fZ(pos.Z()),
    fPath(path)
{
  ;
}

GFPointPath::GFPointPath(double posX, double posY, double posZ,
                         double path)
: fX(posX),  fY(posY),  fZ(posZ),
  fPath(path)
{
;
}


double
GFPointPath::getDist(const GFPointPath& to) const {
  return sqrt(pow(fX - to.X(), 2) +
              pow(fY - to.Y(), 2) +
              pow(fZ - to.Z(), 2));
}


void GFPointPath::Print() const {
  std::cout << "  GFPointPath at "; TVector3(fX, fY, fZ).Print();
  std::cout << "   Path to next point = "<< fPath <<" cm \n";
}

