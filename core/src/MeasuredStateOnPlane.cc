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

#include "TDecompChol.h"

namespace genfit {

void MeasuredStateOnPlane::Print(Option_t*) const {
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

void MeasuredStateOnPlane::blowUpCov(double blowUpFac, bool resetOffDiagonals, double maxVal) {

  unsigned int dim = cov_.GetNcols();

  if (resetOffDiagonals) {
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

  // limit
  if (maxVal > 0.)
    for (unsigned int i=0; i<dim; ++i) {
      for (unsigned int j=0; j<dim; ++j) {
        cov_(i,j) = std::min(cov_(i,j), maxVal);
      }
    }

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
  // For ease of understanding, here's a very explicit implementation
  // that uses the textbook algorithm:
  TMatrixDSym fCovInv, bCovInv, smoothed_cov;
  tools::invertMatrix(forwardState.getCov(), fCovInv);
  tools::invertMatrix(backwardState.getCov(), bCovInv);

  tools::invertMatrix(fCovInv + bCovInv, smoothed_cov);  // one temporary TMatrixDSym

  MeasuredStateOnPlane retVal(forwardState);
  retVal.setState(smoothed_cov*(fCovInv*forwardState.getState() + bCovInv*backwardState.getState())); // four temporary TVectorD's
  retVal.setCov(smoothed_cov);
  return retVal;
#endif

  // This is a numerically stable implementation of the averaging
  // process.  We write S1, S2 for the upper diagonal square roots
  // (Cholesky decompositions) of the covariance matrices, such that
  // C1 = S1' S1 (transposition indicated by ').
  //
  // Then we can write
  //  (C1^-1 + C2^-1)^-1 = (S1inv' S1inv + S2inv' S2inv)^-1
  //                     = ( (S1inv', S2inv') . ( S1inv ) )^-1
  //                       (                    ( S2inv ) )
  //                     = ( (R', 0) . Q . Q' . ( R ) )^-1
  //                       (                    ( 0 ) )
  // where Q is an orthogonal matrix chosen such that R is upper diagonal.
  // Since Q'.Q = 1, this reduces to
  //                     = ( R'.R )^-1
  //                     = R^-1 . (R')^-1.
  // This gives the covariance matrix of the average and its Cholesky
  // decomposition.
  //
  // In order to get the averaged state (writing x1 and x2 for the
  // states) we proceed from
  //  C1^-1.x1 + C2^-1.x2 = (S1inv', S2inv') . ( S1inv.x1 )
  //                                           ( S2inv.x2 )
  // which by the above can be written as
  //                      =  (R', 0) . Q . ( S1inv.x1 )
  //                                       ( S2inv.x2 )
  // with the same (R, Q) as above.
  //
  // The average is then after obvious simplification
  //   average = R^-1 . Q . (S1inv.x1)
  //                        (S2inv.x2)
  //
  // This is what's implemented below, where we make use of the
  // tridiagonal shapes of the various matrices when multiplying or
  // inverting.
  //
  // This turns out not only more numerically stable, but because the
  // matrix operations are simpler, it is also faster than the
  // straightoforward implementation.
  //
  // This is an application of the technique of Golub, G.,
  // Num. Math. 7, 206 (1965) to the least-squares problem underlying
  // averaging.
  TDecompChol d1(forwardState.getCov());
  d1.Decompose();
  TDecompChol d2(backwardState.getCov());
  d2.Decompose();

  int nRows = d1.GetU().GetNrows();
  assert(nRows == d2.GetU().GetNrows());
  TMatrixD S1inv, S2inv;
  tools::transposedInvert(d1.GetU(), S1inv);
  tools::transposedInvert(d2.GetU(), S2inv);

  TMatrixD A(2*nRows, nRows);
  TVectorD b(2 * nRows);
  double *const bk = b.GetMatrixArray();
  double *const Ak = A.GetMatrixArray();
  const double* S1invk = S1inv.GetMatrixArray();
  const double* S2invk = S2inv.GetMatrixArray();
  // S1inv and S2inv are lower triangular.
  for (int i = 0; i < nRows; ++i) {
    double sum1 = 0;
    double sum2 = 0;
    for (int j = 0; j <= i; ++j) {
      Ak[i*nRows + j] = S1invk[i*nRows + j];
      Ak[(i + nRows)*nRows + j] = S2invk[i*nRows + j];
      sum1 += S1invk[i*nRows + j]*forwardState.getState().GetMatrixArray()[j];
      sum2 += S2invk[i*nRows + j]*backwardState.getState().GetMatrixArray()[j];
    }
    bk[i] = sum1;
    bk[i + nRows] = sum2;
  }

  tools::QR(A, b);
  A.ResizeTo(nRows, nRows);

  TMatrixD inv;
  tools::transposedInvert(A, inv);
  b.ResizeTo(nRows);
  for (int i = 0; i < nRows; ++i) {
    double sum = 0;
    for (int j = i; j < nRows; ++j) {
      sum += inv.GetMatrixArray()[j*nRows+i] * b[j];
    }
    b.GetMatrixArray()[i] = sum;
  }
  return MeasuredStateOnPlane(b,
			      TMatrixDSym(TMatrixDSym::kAtA, inv),
			      forwardState.getPlane(),
			      forwardState.getRep(),
			      forwardState.getAuxInfo());
}


} /* End of namespace genfit */
