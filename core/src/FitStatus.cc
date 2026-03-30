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

  value |= opt.Contains("C") ? static_cast<int>(EFields::C) : 0;
  value |= opt.Contains("F") ? static_cast<int>(EFields::F) : 0;
  value |= opt.Contains("L") ? static_cast<int>(EFields::L) : 0;
  value |= opt.Contains("W") ? static_cast<int>(EFields::W) : 0;
  value |= opt.Contains("R") ? static_cast<int>(EFields::R) : 0;
  value |= opt.Contains("M") ? static_cast<int>(EFields::M) : 0;
  value |= opt.Contains("I") ? static_cast<int>(EFields::I) : 0;
  value |= opt.Contains("U") ? static_cast<int>(EFields::U) : 0;
}


bool PruneFlags::hasFlags(Option_t* option) const {
  TString opt = option;
  opt.ToUpper();

  return !((!(value & static_cast<int>(EFields::C)) && opt.Contains("C"))
	   || (!(value & static_cast<int>(EFields::F)) && opt.Contains("F"))
	   || (!(value & static_cast<int>(EFields::L)) && opt.Contains("L"))
	   || (!(value & static_cast<int>(EFields::W)) && opt.Contains("W"))
	   || (!(value & static_cast<int>(EFields::R)) && opt.Contains("R"))
	   || (!(value & static_cast<int>(EFields::M)) && opt.Contains("M"))
	   || (!(value & static_cast<int>(EFields::I)) && opt.Contains("I"))
	   || (!(value & static_cast<int>(EFields::U)) && opt.Contains("U")));
}


bool PruneFlags::isPruned() const {
  return !!value;
}


void PruneFlags::Print(const Option_t*) const {
  printOut << "PruneFlags: ";
  if (value & static_cast<int>(EFields::C)) printOut << "C";
  if (value & static_cast<int>(EFields::F)) printOut << "F";
  if (value & static_cast<int>(EFields::L)) printOut << "L";
  if (value & static_cast<int>(EFields::W)) printOut << "W";
  if (value & static_cast<int>(EFields::R)) printOut << "R";
  if (value & static_cast<int>(EFields::M)) printOut << "M";
  if (value & static_cast<int>(EFields::I)) printOut << "I";
  if (value & static_cast<int>(EFields::U)) printOut << "U";
  printOut << "\n";
}



void FitStatus::Print(const Option_t*) const
{
  printOut << "fitStatus \n";
  if (isFitted_) {
    printOut << " track has been fitted,";
    if (isFitConvergedFully_) {
      printOut << " fit has converged fully,";
    } else if (isFitConvergedPartially_) {
      printOut << " fit has converged partially,";
    } else {
      printOut << " fit has NOT converged,";
    }
    printOut << " " << nFailedPoints_ << " TrackPoints could not be processed,";
    if (trackHasChanged_) {
      printOut << " track has changed since the fit,";
    }
    printOut << " fitted charge = " << charge_ << ", ";
    pruneFlags_.Print();
  }
  else {
    printOut << " track has NOT been fitted,";
  }
}

} /* End of namespace genfit */
