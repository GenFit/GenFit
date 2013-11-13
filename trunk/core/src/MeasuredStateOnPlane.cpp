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

#include "MeasuredStateOnPlane.h"
#include "AbsTrackRep.h"
#include "Exception.h"
#include "Tools.h"

#include <cassert>
#include <iostream>

namespace genfit {

void MeasuredStateOnPlane::Print(Option_t* option) const {
  std::cout << "genfit::MeasuredStateOnPlane ";
  std::cout << "my address " << this << " my plane's address " << this->sharedPlane_.get() << "; use count: " << sharedPlane_.use_count() << std::endl;
  std::cout << " state vector: "; state_.Print();
  std::cout << " covariance matrix: "; cov_.Print();
  if (sharedPlane_ != NULL) {
    std::cout << " defined in plane "; sharedPlane_->Print();
    TVector3 pos, mom;
    TMatrixDSym cov(6,6);
    getRep()->getPosMomCov(*this, pos, mom, cov);
    std::cout << " 3D position: "; pos.Print();
    std::cout << " 3D momentum: "; mom.Print();
    //std::cout << " 6D covariance: "; cov.Print();
  }
}

void MeasuredStateOnPlane::blowUpCov(double blowUpFac, bool resetOffDiagonals) {

  if (resetOffDiagonals) {
    unsigned int dim = cov_.GetNcols();
    for (unsigned int i=0; i<dim; ++i) {
      for (unsigned int j=0; j<dim; ++j) {
        if (i != j)
          cov_(i,j) = 0; // reset off-diagonals
        else
          cov_(i,j) *= blowUpFac; // blow up diagonals
      }
    }
  }
  else
    cov_ *= blowUpFac;

}


MeasuredStateOnPlane calcAverageState(const MeasuredStateOnPlane& forwardState, const MeasuredStateOnPlane& backwardState) {
  // check if both states are defined in the same plane
  if (forwardState.getPlane() != backwardState.getPlane()) {
    Exception e("KalmanFitterInfo::calcAverageState: forwardState and backwardState are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }

  // This code is called a lot, so some effort has gone into reducing
  // the number of temporary objects being constructed.

#if 0
  // For ease of understanding, here's a version of the code that is a
  // few percent slower, depending on the exact use-case, but which
  // makes the math more explicit:
  TMatrixDSym fCovInv, bCovInv, smoothed_cov;
  tools::invertMatrix(forwardState.getCov(), fCovInv);
  tools::invertMatrix(backwardState.getCov(), bCovInv);

  tools::invertMatrix(fCovInv + bCovInv, smoothed_cov);  // one temporary TMatrixDSym

  MeasuredStateOnPlane retVal(forwardState);
  retVal.setState(smoothed_cov*(fCovInv*forwardState.getState() + bCovInv*backwardState.getState())); // four temporary TVectorD's
  retVal.setCov(smoothed_cov);
  return retVal;
#endif

  static TMatrixDSym fCovInv, bCovInv;  // Static to avoid re-constructing for every call
  tools::invertMatrix(forwardState.getCov(), fCovInv);
  tools::invertMatrix(backwardState.getCov(), bCovInv);

  // Using a StateOnPlane here turned out at least as fast as
  // constructing a MeasuredStateOnPlane here, and then resetting its
  // covariance matrix below, even though it means another copy in the
  // return statement.
  StateOnPlane sop(forwardState); // copies auxInfo, plane, rep in the process
                                  // Using 'static' + subsequent assignment is measurably slower.
                                  // Surprisingly, using only TVectorD
                                  // sop(forwardState.getState()) with according
                                  // changes below measured slower.

  sop.getState() *= fCovInv;
  fCovInv += bCovInv;
  tools::invertMatrix(fCovInv);  // This is now the covariance of the average.
  sop.getState() += bCovInv*backwardState.getState();  // one temporary TVectorD
  sop.getState() *= fCovInv;

  return MeasuredStateOnPlane(sop, fCovInv);
}


} /* End of namespace genfit */
