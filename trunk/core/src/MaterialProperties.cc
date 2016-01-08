/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "MaterialProperties.h"

#include "IO.h"

namespace genfit {

bool operator== (const MaterialProperties& lhs, const MaterialProperties& rhs){
  if (&lhs == &rhs)
    return true;
  if (lhs.density_ != rhs.density_ ||
      lhs.Z_ != rhs.Z_ ||
      lhs.A_ != rhs.A_ ||
      lhs.radiationLength_ != rhs.radiationLength_ ||
      lhs.mEE_ != rhs.mEE_)
    return false;

  return true;
}

bool operator!= (const MaterialProperties& lhs, const MaterialProperties& rhs) {
  return !(lhs==rhs);
}


void MaterialProperties::getMaterialProperties(double& density,
                                               double& Z,
                                               double& A,
                                               double& radiationLength,
                                               double& mEE) const {
  density = density_;
  Z = Z_;
  A = A_;
  radiationLength = radiationLength_;
  mEE = mEE_;
}


void MaterialProperties::setMaterialProperties(const double& density,
                                               const double& Z,
                                               const double& A,
                                               const double& radiationLength,
                                               const double& mEE) {
  density_ = density;
  Z_ = Z;
  A_ = A;
  radiationLength_ = radiationLength;
  mEE_ = mEE;
}


void MaterialProperties::Print(const Option_t*) const {
  printOut << "Density = " << density_ << ", \t"
           << "Z = " << Z_ << ", \t"
           << "A = " << A_ << ", \t"
           << "radiationLength = " << radiationLength_ << ", \t"
           << "mEE = " << mEE_ << "\n";
}


} /* End of namespace genfit */
