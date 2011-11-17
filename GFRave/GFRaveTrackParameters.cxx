/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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
#include "GFException.h"

#include <iostream>

using namespace std;

GFRaveTrackParameters::GFRaveTrackParameters() :
  fOriginalTrack(NULL),
  fWeight(0),
  fState(1,6),
  fCov(6,6)
{
  //std::cerr << "GFRaveTrackParameters::GFRaveTrackParameters() => calling default constructor\n";
}


GFRaveTrackParameters::GFRaveTrackParameters(GFTrack* track, double weight, TMatrixT<double> state6, TMatrixT<double> cov6x6) :
  fOriginalTrack(track),
  fWeight(weight),
  fState(state6),
  fCov(cov6x6)
{
  if (fState.GetNrows()!=1 || fState.GetNcols()!=6) {
    GFException exc("GFRaveTrackParameters ==> State is not 1x6!",__LINE__,__FILE__);
    throw exc;
  }
  if (fCov.GetNrows()!=6 || fCov.GetNcols()!=6) {
    GFException exc("GFRaveTrackParameters ==> Covariance is not 6x6!",__LINE__,__FILE__);
    throw exc;
  }
}


TVector3
GFRaveTrackParameters::getPos() const {
  return TVector3(fState[0][0], fState[0][1], fState[0][2]);
}


TVector3
GFRaveTrackParameters::getMom() const {
  return TVector3(fState[0][3], fState[0][4], fState[0][5]);
}


double
GFRaveTrackParameters::getCharge() const {
  return fOriginalTrack->getCardinalRep()->getCharge();
}


double
GFRaveTrackParameters::getPdg() const{
  return fOriginalTrack->getCardinalRep()->getPDG();
}


void
GFRaveTrackParameters::Print(const Option_t*) const {
  std::cout << "weight: " << getWeight() << "\n";
  std::cout << "state: "; getState().Print();
  std::cout << "cov: "; getCov().Print();
  std::cout << "GFTrack: "; getTrack()->Print();
}
