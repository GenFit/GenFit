/* Copyright 2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Johannes Rauch & Tobias Schlüter

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


#include "FitStatus.h"

#include <iostream>

namespace genfit {

void FitStatus::Print(const Option_t*) const
{
  std::cout << "fitStatus \n";
  if (isFitted_) {
    std::cout << " track has been fitted,";
    if (isFitConverged_)
      std::cout << " fit has converged,";
    else
      std::cout << " fit has NOT converged,";
    if (trackHasChanged_) std::cout << " track has changed since the fit,";
    if (trackIsPruned_) std::cout << " track is pruned,";
    std::cout << " fitted charge = " << charge_ << " \n";
  }
  else
    std::cout << " track has NOT been fitted,";
}

} /* End of namespace genfit */
