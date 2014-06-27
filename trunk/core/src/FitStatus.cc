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

  C = opt.Contains("C") || C;
  F = opt.Contains("F") || F;
  L = opt.Contains("L") || L;
  W = opt.Contains("W") || W;
  R = opt.Contains("R") || R;
  M = opt.Contains("M") || M;
  I = opt.Contains("I") || I;
  U = opt.Contains("U") || U;
}


bool PruneFlags::hasFlags(Option_t* option) const {
  TString opt = option;
  opt.ToUpper();

  if ((opt.Contains("C") && !C) ||
      (opt.Contains("F") && !F) ||
      (opt.Contains("L") && !L) ||
      (opt.Contains("W") && !W) ||
      (opt.Contains("R") && !R) ||
      (opt.Contains("M") && !M) ||
      (opt.Contains("I") && !I) ||
      (opt.Contains("U") && !U) )
    return false;

  return true;
}


bool PruneFlags::isPruned() const {
  return (C || F || L || W || R || M || I || U);
}


void PruneFlags::Print(const Option_t* option) const {
  std::cout << "PruneFlags: ";
  if (C) std::cout << "C";
  if (F) std::cout << "F";
  if (L) std::cout << "L";
  if (W) std::cout << "W";
  if (R) std::cout << "R";
  if (M) std::cout << "M";
  if (I) std::cout << "I";
  if (U) std::cout << "U";
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
    std::cout << " fitted charge = " << charge_;
    pruneFlags_.Print();
  }
  else
    std::cout << " track has NOT been fitted,";
}

} /* End of namespace genfit */
