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


#include "KalmanFitStatus.h"

#include "IO.h"

namespace genfit {

void KalmanFitStatus::Print(const Option_t*) const
{
  FitStatus::Print();
  if (fittedWithDaf_) printOut << " track has been fitted with DAF,";
  if (fittedWithReferenceTrack_) printOut << " track has been fitted with reference track,";
  if (isFitted_) {
    printOut << " numIterations = " << numIterations_ << ", ";
    printOut << "track length = " << trackLen_ << ", ";
    printOut << "fChi2 = " << fChi2_ << ", ";
    printOut << "bChi2 = " << FitStatus::getChi2() << ", ";
    printOut << "fNdf = " << fNdf_ << ", ";
    printOut << "bNdf = " << FitStatus::getNdf() << ", ";
    printOut << "fPVal = " << getForwardPVal() << ", ";
    printOut << "bPVal = " << getBackwardPVal() << "\n";
  }
  printOut << "\n";
}

} /* End of namespace genfit */
