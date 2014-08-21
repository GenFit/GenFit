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

#include <TString.h>

#include <iostream>

namespace genfit {

PruneFlags::PruneFlags() {
  reset();
}


void PruneFlags::reset() {
  memset(this, 0, sizeof *this);
}


void PruneFlags::setFlags(Option_t* option) {
  TString opt = option;
  opt.ToUpper();

  value |= opt.Contains("C") ? C : 0;
  value |= opt.Contains("F") ? F : 0;
  value |= opt.Contains("L") ? L : 0;
  value |= opt.Contains("W") ? W : 0;
  value |= opt.Contains("R") ? R : 0;
  value |= opt.Contains("M") ? M : 0;
  value |= opt.Contains("I") ? I : 0;
  value |= opt.Contains("U") ? U : 0;
}


bool PruneFlags::hasFlags(Option_t* option) const {
  TString opt = option;
  opt.ToUpper();

  return !((!(value & C) && opt.Contains("C"))
	   || (!(value & F) && opt.Contains("F"))
	   || (!(value & L) && opt.Contains("L"))
	   || (!(value & W) && opt.Contains("W"))
	   || (!(value & R) && opt.Contains("R"))
	   || (!(value & M) && opt.Contains("M"))
	   || (!(value & I) && opt.Contains("I"))
	   || (!(value & U) && opt.Contains("U")));
}


bool PruneFlags::isPruned() const {
  return !!value;
}


void PruneFlags::Print(const Option_t*) const {
  std::cout << "PruneFlags: ";
  if (value & C) std::cout << "C";
  if (value & F) std::cout << "F";
  if (value & L) std::cout << "L";
  if (value & W) std::cout << "W";
  if (value & R) std::cout << "R";
  if (value & M) std::cout << "M";
  if (value & I) std::cout << "I";
  if (value & U) std::cout << "U";
  std::cout << "\n";
}



void FitStatus::Print(const Option_t*) const
{
  std::cout << "fitStatus \n";
  if (isFitted_) {
    std::cout << " track has been fitted,";
    if (isFitConvergedFully_)
      std::cout << " fit has converged fully,";
    else if (isFitConvergedPartially_)
      std::cout << " fit has converged partially,";
    else
      std::cout << " fit has NOT converged,";
    std::cout << " " << nFailedPoints_ << " TrackPoints could not be processed,";
    if (trackHasChanged_) std::cout << " track has changed since the fit,";
    std::cout << " fitted charge = " << charge_ << ", ";
    pruneFlags_.Print();
  }
  else
    std::cout << " track has NOT been fitted,";
}

} /* End of namespace genfit */
