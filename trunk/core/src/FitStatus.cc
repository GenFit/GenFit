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
#include "IO.h"

#include <TString.h>

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
  printOut << "PruneFlags: ";
  if (value & C) printOut << "C";
  if (value & F) printOut << "F";
  if (value & L) printOut << "L";
  if (value & W) printOut << "W";
  if (value & R) printOut << "R";
  if (value & M) printOut << "M";
  if (value & I) printOut << "I";
  if (value & U) printOut << "U";
  printOut << "\n";
}



void FitStatus::Print(const Option_t*) const
{
  printOut << "fitStatus \n";
  if (isFitted_) {
    printOut << " track has been fitted,";
    if (isFitConvergedFully_)
      printOut << " fit has converged fully,";
    else if (isFitConvergedPartially_)
      printOut << " fit has converged partially,";
    else
      printOut << " fit has NOT converged,";
    printOut << " " << nFailedPoints_ << " TrackPoints could not be processed,";
    if (trackHasChanged_) printOut << " track has changed since the fit,";
    printOut << " fitted charge = " << charge_ << ", ";
    pruneFlags_.Print();
  }
  else
    printOut << " track has NOT been fitted,";
}

} /* End of namespace genfit */
