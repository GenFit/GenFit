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
  fCharge(0),
  fPdg(0)
{

}


GFRaveTrackParameters::GFRaveTrackParameters(TMatrixT<double> state, TMatrixT<double> cov, double charge, int pdg) :
  fState(state),
  fCov(cov),
  fCharge(charge),
  fPdg(pdg)
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
