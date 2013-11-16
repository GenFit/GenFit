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


#include "GFRaveMagneticField.h"
#include <FieldManager.h>

#include <iostream>


namespace genfit {

GFRaveMagneticField *
GFRaveMagneticField::copy() const{
  return new GFRaveMagneticField(*this);
}


rave::Vector3D
GFRaveMagneticField::inTesla ( const rave::Point3D & position) const
{
  TVector3 pos(position.x(), position.y(), position.z());

  TVector3 B = FieldManager::getInstance()->getFieldVal(pos); // magnetic field in kGauss
  B *= 1.E-1;

  return rave::Vector3D (B.X(), B.Y(), B.Z());
}


} /* End of namespace genfit */
