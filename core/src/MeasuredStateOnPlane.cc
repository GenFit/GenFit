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
#include "IO.h"
#include "RootEigenTransformations.h"
#include "EigenMatrixTypedefs.h"

#include <cassert>

#include "TDecompChol.h"

namespace genfit {

void MeasuredStateOnPlane::Print(Option_t*) const {
  printOut << "genfit::MeasuredStateOnPlane ";
  printOut << "my address " << this << " my plane's address " << this->sharedPlane_.get() << "; use count: " << sharedPlane_.use_count() << std::endl;
  printOut << " state vector: "; state_.Print();
  printOut << " covariance matrix: "; cov_.Print();
  if (sharedPlane_ != nullptr) {
    printOut << " defined in plane "; sharedPlane_->Print();
    TVector3 pos, mom;
    TMatrixDSym cov(6,6);
    getRep()->getPosMomCov(*this, pos, mom, cov);
    printOut << " 3D position: "; pos.Print();
    printOut << " 3D momentum: "; mom.Print();
    //printOut << " 6D covariance: "; cov.Print();
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


MeasuredStateOnPlane calcAverageState(const MeasuredStateOnPlane &forwardState, const MeasuredStateOnPlane &backwardState) {
    // check if both states are defined in the same plane
    if (forwardState.getPlane() != backwardState.getPlane()) {
        Exception e(
                "KalmanFitterInfo::calcAverageState: forwardState and backwardState are not defined in the same plane.",
                __LINE__, __FILE__);
        throw e;
    }

//#define TEXTBOOK
#ifdef TEXTBOOK
    // For ease of understanding, here's a very explicit implementation
    // that uses the textbook algorithm:
    const Vector5 fState(rootVectorToEigenVector<5>(forwardState.getState()));
    const Vector5 bState(rootVectorToEigenVector<5>(backwardState.getState()));
    const Matrix5x5Sym fCovInv(rootMatrixSymToEigenMatrix<5>(forwardState.getCov()).inverse());
    const Matrix5x5Sym bCovInv(rootMatrixSymToEigenMatrix<5>(backwardState.getCov()).inverse());
    const Matrix5x5Sym smoothedCov((fCovInv + bCovInv).inverse());

    return MeasuredStateOnPlane(
            eigenVectorToRootVector<5>(smoothedCov*(fCovInv * fState + bCovInv * bState)),
            eigenMatrixToRootMatrixSym<5>(smoothedCov),
            forwardState.getPlane(),
            forwardState.getRep(),
            forwardState.getAuxInfo()
    );
#else

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
    Eigen::LLT<Eigen::MatrixXd> lltFwState(rootMatrixSymToEigenMatrix<5>(forwardState.getCov()));
    if (lltFwState.info() == Eigen::NumericalIssue) {
        throw Exception("KalmanFitterInfo::calcAverageState: ill-conditioned covariance matrix of forward state.",
                        __LINE__, __FILE__);
    }
    Eigen::LLT<Eigen::MatrixXd> lltBwState(rootMatrixSymToEigenMatrix<5>(backwardState.getCov()));
    if (lltFwState.info() == Eigen::NumericalIssue) {
        throw Exception("KalmanFitterInfo::calcAverageState: ill-conditioned covariance matrix of backward state.",
                        __LINE__, __FILE__);
    }
    const Matrix5x5 S1inv(lltFwState.matrixU().transpose().solve(Matrix5x5::Identity()));
    const Matrix5x5 S2inv(lltBwState.matrixU().transpose().solve(Matrix5x5::Identity()));
    Eigen::Matrix<Scalar, 10, 5> A;
    A.block<5, 5>(0, 0) = S1inv;
    A.block<5, 5>(5, 0) = S2inv;
    Eigen::Matrix<Scalar, 10, 1> b;
    b.block<5, 1>(0, 0) = S1inv * rootVectorToEigenVector<5>(forwardState.getState());
    b.block<5, 1>(5, 0) = S2inv * rootVectorToEigenVector<5>(backwardState.getState());

    Eigen::ColPivHouseholderQR<decltype(A)> QRdecomp(A);
    const Matrix5x5 R = QRdecomp.matrixR().block<5, 5>(0, 0);
    const Vector5 result = QRdecomp.solve(b).block<5, 1>(0, 0);

    Matrix5x5 inv(R.triangularView<Eigen::Lower>().solve(Matrix5x5::Identity()));

    return MeasuredStateOnPlane(eigenVectorToRootVector<5>(result),
                                eigenMatrixToRootMatrixSym<5>(inv.transpose() * inv),
                                forwardState.getPlane(),
                                forwardState.getRep(),
                                forwardState.getAuxInfo());
#endif
}

} /* End of namespace genfit */