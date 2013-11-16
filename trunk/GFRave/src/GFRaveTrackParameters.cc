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


#include "GFRaveTrackParameters.h"
#include "GFRaveConverters.h"
#include "Exception.h"

#include <iostream>


namespace genfit {


GFRaveTrackParameters::GFRaveTrackParameters() :
  originalTrack_(NULL),
  weight_(0),
  state_(6),
  cov_(6,6),
  hasSmoothedData_(false)
{
  ;
}


GFRaveTrackParameters::GFRaveTrackParameters(const Track* track, MeasuredStateOnPlane* originalState, double weight, const TVectorD & state6, const TMatrixDSym & cov6x6, bool isSmoothed) :
  originalTrack_(const_cast<Track*>(track)),
  weight_(weight),
  state_(state6),
  cov_(cov6x6),
  hasSmoothedData_(isSmoothed)
{
  if (state_.GetNrows() != 6) {
    Exception exc("GFRaveTrackParameters ==> State is not 6D!",__LINE__,__FILE__);
    throw exc;
  }
  if (cov_.GetNrows()!=6) {
    Exception exc("GFRaveTrackParameters ==> Covariance is not 6D!",__LINE__,__FILE__);
    throw exc;
  }

}


GFRaveTrackParameters::GFRaveTrackParameters(const Track* track, MeasuredStateOnPlane* originalState, double weight) :
  originalTrack_(const_cast<Track*>(track)),
  weight_(weight),
  state_(1,6),
  cov_(6,6),
  hasSmoothedData_(false)
{
  ;
}


TVector3
GFRaveTrackParameters::getPos() const {
  return TVector3(state_[0], state_[1], state_[2]);
}


TVector3
GFRaveTrackParameters::getMom() const {
  return TVector3(state_[3], state_[4], state_[5]);
}


double
GFRaveTrackParameters::getCharge() const {
  return getTrack()->getFitStatus()->getCharge();
}


double
GFRaveTrackParameters::getPdg() const{
  if (hasTrack())
    return getTrack()->getCardinalRep()->getPDG();
  else {
    Exception exc("GFRaveTrackParameters::getPdg() ==> no genfit::Track available!",__LINE__,__FILE__);
    throw exc;
  }
}


void
GFRaveTrackParameters::Print(const Option_t*) const {
  std::cout << "weight: " << getWeight() << "\n";
  if (!hasSmoothedData_) std::cout << "state and cov are NOT smoothed! \n";
  std::cout << "state: "; getState().Print();
  std::cout << "cov: "; getCov().Print();
  if (hasTrack()) {std::cout << "genfit::Track: "; getTrack()->Print();}
  else std::cout << "NO genfit::Track pointer \n";
}


} /* End of namespace genfit */
