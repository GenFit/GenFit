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

#include "StepLimits.h"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <limits>


namespace genfit {

const double StepLimits::maxLimit_ = 99.E99;


StepLimits& StepLimits::operator=(const StepLimits& other) {
  for (unsigned int i=1; i<ENUM_NR_ITEMS; ++i) {
    limits_[i] = other.limits_[i];
  }

  stepSign_ = other.stepSign_;

  return *this;
}


std::pair<StepLimitType, double> StepLimits::getLowestLimit(double margin) const {

  double lowest(maxLimit_);
  unsigned int iLowest(0);

  for (unsigned int i=1; i<ENUM_NR_ITEMS; ++i) {

    // lowest hard limit may exceed lowest soft limit by up to #margin
    if (i == int(stp_sMaxArg))
      lowest *= (1+margin);

    if (limits_[i] < lowest) {
      lowest = limits_[i];
      iLowest = i;
    }
  }

  return std::pair<StepLimitType, double>(static_cast<StepLimitType>(iLowest), lowest);
}


double StepLimits::getLowestLimitVal(double margin) const {

  double lowest(maxLimit_);

  for (unsigned int i=1; i<ENUM_NR_ITEMS; ++i) {

    // lowest hard limit may exceed lowest soft limit by up to #margin
    if (i == int(stp_sMaxArg))
      lowest *= (1+margin);

    if (limits_[i] < lowest) {
      lowest = limits_[i];
    }
  }

  return lowest;
}


void StepLimits::reduceLimit(StepLimitType type, double value) {
  assert (type != stp_noLimit);
  value = fabs(value);

  if (limits_[type] > value)
    limits_[type] = value;
}


void StepLimits::setStepSign(char signedVal) {
  if (signedVal < 0)
    stepSign_ = -1;
  else
    stepSign_ = 1;
}

void StepLimits::setStepSign(double signedVal) {
  if (signedVal < 0.)
    stepSign_ = -1;
  else
    stepSign_ = 1;
}


void StepLimits::reset() {
  for (unsigned int i=1; i<ENUM_NR_ITEMS; ++i) {
    limits_[i] = maxLimit_;
  }
  stepSign_ = 1;
}


void StepLimits::Print() {
  for (unsigned int i=0; i<ENUM_NR_ITEMS; ++i) {
    if (limits_[i] >= maxLimit_)
      continue;

    std::cout << "   | " << limits_[i] << " cm due to ";
    switch (static_cast<StepLimitType>(i)) {
    case stp_noLimit:
      break;
    case stp_fieldCurv:
      std::cout << "stp_fieldCurv (medium limit): stepsize limited by curvature and magnetic field inhomogenities";
      break;
    case stp_momLoss:
      std::cout << "stp_momLoss (medium limit): stepsize limited by stepper because maximum momLoss is reached";
      break;
    case stp_sMax:
      std::cout << "stp_sMax (medium limit): stepsize limited by SMax defined in #estimateStep()";
      break;
    case stp_sMaxArg:
      std::cout << "stp_sMaxArg (hard limit): stepsize limited by argument maxStepArg passed to #estimateStep()";
      break;
    case stp_boundary:
      std::cout << "stp_boundary (hard limit): stepsize limited by stepper because material boundary is encountered";
      break;
    case stp_plane:
      std::cout << "stp_plane (hard limit):  stepsize limited because destination plane is reached";
      break;
    case ENUM_NR_ITEMS:
      break;
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

} /* End of namespace genfit */
