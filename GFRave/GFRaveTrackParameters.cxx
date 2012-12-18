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
  fOriginalRep(NULL),
  fWeight(0),
  fState(1,6),
  fCov(6,6),
  fHasSmoothedData(false)
{
  ;
}


GFRaveTrackParameters::GFRaveTrackParameters(GFTrack* track, GFAbsTrackRep* rep, double weight, const TVectorD & state6, const TMatrixDSym & cov6x6, bool isSmoothed) :
  fOriginalTrack(track),
  fOriginalRep(rep),
  fWeight(weight),
  fState(state6),
  fCov(cov6x6),
  fHasSmoothedData(isSmoothed)
{
  if (fState.GetNrows() != 6) {
    GFException exc("GFRaveTrackParameters ==> State is not 6D!",__LINE__,__FILE__);
    throw exc;
  }
  if (fCov.GetNrows()!=6) {
    GFException exc("GFRaveTrackParameters ==> Covariance is not 6D!",__LINE__,__FILE__);
    throw exc;
  }
}


GFRaveTrackParameters::GFRaveTrackParameters(GFTrack* track, GFAbsTrackRep* rep, double weight) :
  fOriginalTrack(track),
  fOriginalRep(rep),
  fWeight(weight),
  fState(1,6),
  fCov(6,6),
  fHasSmoothedData(false)
{
  ;
}


TVector3
GFRaveTrackParameters::getPos() const {
  return TVector3(fState[0], fState[1], fState[2]);
}


TVector3
GFRaveTrackParameters::getMom() const {
  return TVector3(fState[3], fState[4], fState[5]);
}


double
GFRaveTrackParameters::getCharge() const {
  return getTrack()->getCardinalRep()->getCharge();
}


double
GFRaveTrackParameters::getPdg() const{
  return getTrack()->getCardinalRep()->getPDG();
}


void
GFRaveTrackParameters::Print(const Option_t*) const {
  std::cout << "weight: " << getWeight() << "\n";
  if (!fHasSmoothedData) std::cout << "state and cov are NOT smoothed! \n";
  std::cout << "state: "; getState().Print();
  std::cout << "cov: "; getCov().Print();
  if (hasTrack()) {std::cout << "GFTrack: "; getTrack()->Print();}
  else std::cout << "NO GFTrack pointer \n";
  if (hasRep()) {std::cout << "GFAbsTrackRep: "; getRep()->Print();}
  else std::cout << "NO GFAbsTrackRep pointer \n";
}
