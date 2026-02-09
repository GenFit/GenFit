/*
 * GblPoint.cpp
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \file
 *  GblPoint methods.
 *
 *  \author Claus Kleinwort, DESY, 2011 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2011 - 2023 Deutsches Elektronen-Synchroton,
 *  Member of the Helmholtz Association, (DESY), HAMBURG, GERMANY \n\n
 *  This library is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version. \n\n
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details. \n\n
 *  You should have received a copy of the GNU Library General Public
 *  License along with this program (see the file COPYING.LIB for more
 *  details); if not, write to the Free Software Foundation, Inc.,
 *  675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "GblPoint.h"
using namespace Eigen;

//! Namespace for the general broken lines package
namespace gbl {

/// Create a point.
/**
 * Create point on (initial) trajectory. Needs transformation jacobian from previous point.
 * \param [in] aJacobian Transformation jacobian from previous point
 * \param [in] numMeasReserve number of measurements to reserve (space for)
 */
GblPoint::GblPoint(const Matrix5d &aJacobian, unsigned int numMeasReserve) :
		theLabel(0), theOffset(0), theType(0), p2pJacobian(aJacobian), scatDim(
				0) {
	theMeasurements.reserve(numMeasReserve);
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Create a point.
/**
 * Create point on (initial) trajectory. Needs transformation jacobian from previous point.
 * \param [in] aJacobian Transformation jacobian from previous point
 */
GblPoint::GblPoint(const TMatrixD &aJacobian) :
		theLabel(0), theOffset(0), theType(0), scatDim(0) {

	for (unsigned int i = 0; i < 5; ++i) {
		for (unsigned int j = 0; j < 5; ++j) {
			p2pJacobian(i, j) = aJacobian(i, j);
		}
	}
}
#endif

GblPoint::~GblPoint() {
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Add a measurement to a point.
/**
 * Add measurement (in meas. system) with diagonal precision (inverse covariance) matrix.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aProjection Projection from local to measurement system
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (diagonal)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
void GblPoint::addMeasurement(const TMatrixD &aProjection,
		const TVectorD &aResiduals, const TVectorD &aPrecision,
		double minPrecision) {
	theMeasurements.emplace_back(aProjection, aResiduals, aPrecision,
			minPrecision);
}

/// Add a measurement to a point.
/**
 * Add measurement (in meas. system) with arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aProjection Projection from local to measurement system
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (matrix)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
void GblPoint::addMeasurement(const TMatrixD &aProjection,
		const TVectorD &aResiduals, const TMatrixDSym &aPrecision,
		double minPrecision) {
	theMeasurements.emplace_back(aProjection, aResiduals, aPrecision,
			minPrecision);
}

/// Add a measurement to a point.
/**
 * Add measurement in local system with diagonal precision (inverse covariance) matrix.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (diagonal)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
void GblPoint::addMeasurement(const TVectorD &aResiduals,
		const TVectorD &aPrecision, double minPrecision) {
	theMeasurements.emplace_back(aResiduals, aPrecision, minPrecision);
}

/// Add a measurement to a point.
/**
 * Add measurement in local system with arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (matrix)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
void GblPoint::addMeasurement(const TVectorD &aResiduals,
		const TMatrixDSym &aPrecision, double minPrecision) {
	theMeasurements.emplace_back(aResiduals, aPrecision, minPrecision);
}
#endif

/// Check for measurements at a point.
/**
 * Get number of measurement (0 = none).
 * \return measurements size
 */
unsigned int GblPoint::numMeasurements() const {
	return theMeasurements.size();
}

/// Add a thin scatterer to a point.
/**
 * Add scatterer with arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized. Changes local track direction.
 *
 * The precision matrix for the local slopes is defined by the
 * angular scattering error theta_0 and the scalar products c_1, c_2 of the
 * offset directions in the local frame with the track direction:
 *
 *            (1 - c_1*c_1 - c_2*c_2)   |  1 - c_1*c_1     - c_1*c_2  |
 *       P =  ----------------------- * |                             |
 *                theta_0*theta_0       |    - c_1*c_2   1 - c_2*c_2  |
 *
 * One offset is defined at the point.
 * The comparison with the previous and next offsets allows to define a kink.
 *
 * \tparam Precision   Precision matrix or vector (with diagonal)
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (full matrix)
 */
void GblPoint::addScatterer(const Eigen::Vector2d &aResiduals,
		const Eigen::Matrix2d &aPrecision) {
	scatDim = 2;
	// arbitrary precision matrix
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> scatEigen { aPrecision };
	scatTransformation.topLeftCorner<2, 2>() = scatEigen.eigenvectors();
	scatTransformation.topLeftCorner<2, 2>().transposeInPlace();
	scatResiduals.head<2>() = scatTransformation.topLeftCorner<2, 2>()
			* aResiduals;
	scatPrecision.head<2>() = scatEigen.eigenvalues();
}

/// Add a thin scatterer to a point.
/**
 * Add scatterer with diagonal precision (inverse covariance) matrix.
 * Changes local track direction.
 *
 * The precision matrix for the local slopes is defined by the
 * angular scattering error theta_0 and the scalar products c_1, c_2 of the
 * offset directions in the local frame with the track direction:
 *
 *            (1 - c_1*c_1 - c_2*c_2)   |  1 - c_1*c_1     - c_1*c_2  |
 *       P =  ----------------------- * |                             |
 *                theta_0*theta_0       |    - c_1*c_2   1 - c_2*c_2  |
 *
 * For a diagonal precision matrix at least one of the scalar products c_1, c_2
 * must be zero.
 *
 * One offset is defined at the point.
 * The comparison with the previous and next offsets allows to define a kink.
 *
 * \tparam Precision   Precision matrix or vector (with diagonal)
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (vector with diagonal)
 */
void GblPoint::addScatterer(const Eigen::Vector2d &aResiduals,
		const Eigen::Vector2d &aPrecision) {
	scatDim = 2;
	scatResiduals.head<2>() = aResiduals;
	scatPrecision.head<2>() = aPrecision;
	scatTransformation.topLeftCorner<2, 2>().setIdentity();
}

/// Add a thick scatterer to a point.
/**
 * Add scatterer with arbitrary precision (inverse 4x4 covariance) matrix.
 * Will be diagonalized. Local track direction and position change (correlated).
 *
 * Combination of several scatterers extrapolated to point
 * (yielding dense covariance matrix). E.g. sum of thin scatterers (covariance V_i,
 * jacobian J_i for slopes and offsets):
 *
 *       P = [sum(J_i*V_i*J_i^t)]^-1
 *
 * For each thin scatterer the covarinace matrix V_i is defined by the
 * angular scattering error theta_0,i and the scalar products c_1,i, c_2,i
 * of the offset directions in the local frame with the track direction:
 *
 *                                         | 1 - c_2,i^2  c_1,i*c_2,i     0     0 |
 *                                         |                                      |
 *                     theta_0,i^2         | c_1,i*c_2,i  1 - c_1,i^2     0     0 |
 *       V_i = ------------------------- * |                                      |
 *             (1 - c_1,i^2 - c_2,i^2)^2   |      0            0          0     0 |
 *                                         |                                      |
 *                                         |      0            0          0     0 |
 *
 * Two offsets are defined at the point to describe a step (in the trajectory).
 * The comparison with the previous and next offsets allows to define a kink.
 *
 * \tparam Precision   Precision matrix
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (full 4x4 matrix)
 */
void GblPoint::addThickScatterer(const Eigen::Vector4d &aResiduals,
		const Eigen::Matrix4d &aPrecision) {
	scatDim = 4;
	// dense precision matrix
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> scatEigen { aPrecision };
	scatTransformation = scatEigen.eigenvectors();
	scatTransformation.transposeInPlace();
	scatResiduals = scatTransformation * aResiduals;
	scatPrecision = scatEigen.eigenvalues();
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Add a thin scatterer to a point.
/**
 * Add scatterer with diagonal precision (inverse 2x2 covariance) matrix.
 * Changes local track direction.
 *
 * The precision matrix for the local slopes is defined by the
 * angular scattering error theta_0 and the scalar products c_1, c_2 of the
 * offset directions in the local frame with the track direction:
 *
 *            (1 - c_1*c_1 - c_2*c_2)   |  1 - c_1*c_1     - c_1*c_2  |
 *       P =  ----------------------- * |                             |
 *                theta_0*theta_0       |    - c_1*c_2   1 - c_2*c_2  |
 *
 * For a diagonal precision matrix at least one of the scalar products c_1, c_2
 * must be zero.
 *
 * One offset is defined at the point.
 * The comparison with the previous and next offsets allows to define a kink.
 *
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (diagonal of inverse covariance matrix)
 */
void GblPoint::addScatterer(const TVectorD &aResiduals,
		const TVectorD &aPrecision) {
	scatDim = 2;
	scatResiduals(0) = aResiduals[0];
	scatResiduals(1) = aResiduals[1];
	scatPrecision(0) = aPrecision[0];
	scatPrecision(1) = aPrecision[1];
	scatTransformation.setIdentity();
}

/// Add a thin scatterer to a point.
/**
 * Add scatterer with arbitrary precision (inverse 2x2 covariance) matrix.
 * Will be diagonalized. Changes local track direction.
 *
 * The precision matrix for the local slopes is defined by the
 * angular scattering error theta_0 and the scalar products c_1, c_2 of the
 * offset directions in the local frame with the track direction:
 *
 *            (1 - c_1*c_1 - c_2*c_2)   |  1 - c_1*c_1     - c_1*c_2  |
 *       P =  ----------------------- * |                             |
 *                theta_0*theta_0       |    - c_1*c_2   1 - c_2*c_2  |
 *
 * One offset is defined at the point.
 * The comparison with the previous and next offsets allows to define a kink.
 *
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (matrix)
 */
void GblPoint::addScatterer(const TVectorD &aResiduals,
		const TMatrixDSym &aPrecision) {
	scatDim = 2;
	TMatrixDSymEigen scatEigen(aPrecision);
	TMatrixD aTransformation = scatEigen.GetEigenVectors();
	aTransformation.T();
	TVectorD transResiduals = aTransformation * aResiduals;
	TVectorD transPrecision = scatEigen.GetEigenValues();
	for (unsigned int i = 0; i < 2; ++i) {
		scatResiduals(i) = transResiduals[i];
		scatPrecision(i) = transPrecision[i];
		for (unsigned int j = 0; j < 2; ++j) {
			scatTransformation(i, j) = aTransformation(i, j);
		}
	}
}

/// Add a thick scatterer to a point.
/**
 * Add scatterer with arbitrary precision (inverse 4x4 covariance) matrix.
 * Will be diagonalized. Local track direction and position change (correlated).
 *
 * Combination of several scatterers extrapolated to point
 * (yielding dense covariance matrix). E.g. sum of thin scatterers (covariance V_i,
 * jacobian J_i for slopes and offsets):
 *
 *       P = [sum(J_i*V_i*J_i^t)]^-1
 *
 * For each thin scatterer the covarinace matrix V_i is defined by the
 * angular scattering error theta_0,i and the scalar products c_1,i, c_2,i
 * of the offset directions in the local frame with the track direction:
 *
 *                                         | 1 - c_2,i^2  c_1,i*c_2,i     0     0 |
 *                                         |                                      |
 *                     theta_0,i^2         | c_1,i*c_2,i  1 - c_1,i^2     0     0 |
 *       V_i = ------------------------- * |                                      |
 *             (1 - c_1,i^2 - c_2,i^2)^2   |      0            0          0     0 |
 *                                         |                                      |
 *                                         |      0            0          0     0 |
 *
 * Two offsets are defined at the point to describe a step (in the trajectory).
 * The comparison with the previous and next offsets allows to define a kink.
 *
 * \tparam Precision   Precision matrix
 * \param [in] aResiduals Scatterer residuals
 * \param [in] aPrecision Scatterer precision (full 4x4 matrix)
 */
void GblPoint::addThickScatterer(const TVectorD &aResiduals,
		const TMatrixDSym &aPrecision) {
	scatDim = 4;
	TMatrixDSymEigen scatEigen(aPrecision);
	TMatrixD aTransformation = scatEigen.GetEigenVectors();
	aTransformation.T();
	TVectorD transResiduals = aTransformation * aResiduals;
	TVectorD transPrecision = scatEigen.GetEigenValues();
	for (unsigned int i = 0; i < 4; ++i) {
		scatResiduals(i) = transResiduals[i];
		scatPrecision(i) = transPrecision[i];
		for (unsigned int j = 0; j < 4; ++j) {
			scatTransformation(i, j) = aTransformation(i, j);
		}
	}
}
#endif

/// Get scatterer dimension.
/**
 * Get dimension of scatterer (0 = none, 2 = thin, 4 = thick).
 * \return scatterer dimension
 */
unsigned int GblPoint::getScatDim() const {
	return scatDim;
}

/// Retrieve scatterer of a point.
/**
 * \param [out] aTransformation Scatterer transformation from diagonalization
 * \param [out] aResiduals Scatterer residuals
 * \param [out] aPrecision Scatterer precision (diagonal)
 */
void GblPoint::getScatterer(Matrix4d &aTransformation, Vector4d &aResiduals,
		Vector4d &aPrecision) const {
	aTransformation.topLeftCorner(scatDim, scatDim) =
			scatTransformation.topLeftCorner(scatDim, scatDim);
	aResiduals.head(scatDim) = scatResiduals.head(scatDim);
	aPrecision.head(scatDim) = scatPrecision.head(scatDim);
}

/// Retrieve (reduced) scatterer of a point.
/**
 * Reduce to positional scatterer.
 *
 * \param [out] aTransformation Scatterer transformation from diagonalization
 * \param [out] aResiduals Scatterer residuals
 * \param [out] aPrecision Scatterer precision (diagonal)
 */
void GblPoint::getReducedScatterer(Matrix4d &aTransformation,
		Vector4d &aResiduals, Vector4d &aPrecision) const {
	// get scatterer input back from diagonalization
	Matrix4d scatPrec = scatTransformation.transpose()
			* scatPrecision.asDiagonal() * scatTransformation;
	Vector4d scatRes = scatTransformation.transpose() * scatResiduals;
	// block matrix algebra to get positional precision
	Matrix2d posPrec = scatPrec.bottomRightCorner<2, 2>()
			- scatPrec.bottomLeftCorner<2, 2>()
					* scatPrec.topLeftCorner<2, 2>().inverse()
					* scatPrec.topRightCorner<2, 2>();
	// scatPrec with only positional precision
	scatPrec.setZero();
	scatPrec.bottomRightCorner<2, 2>() = posPrec;
	// diagonalize again
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> scatEigen { scatPrec };
	//
	aTransformation = scatEigen.eigenvectors().transpose();
	aResiduals = aTransformation * scatRes;
	aPrecision = scatEigen.eigenvalues();
}

/// Get scatterer transformation (from diagonalization).
/**
 * \param [out] aTransformation Transformation matrix
 */
void GblPoint::getScatTransformation(MatrixXd &aTransformation) const {
	aTransformation.resize(scatDim, scatDim);
	if (scatDim > 0) {
		aTransformation = scatTransformation.topLeftCorner(scatDim, scatDim);
	} else {
		aTransformation.setIdentity();
	}
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Add local derivatives to a point.
/**
 * Point needs to have a measurement.
 * \param [in] aDerivatives Local derivatives (matrix)
 */
void GblPoint::addLocals(const TMatrixD &aDerivatives) {
	if (theMeasurements.size())
		theMeasurements.back().addLocals(aDerivatives);
}

/// Add global derivatives to a point.
/**
 * Point needs to have a measurement.
 * \param [in] aLabels Global derivatives labels
 * \param [in] aDerivatives Global derivatives (matrix)
 */
void GblPoint::addGlobals(const std::vector<int> &aLabels,
		const TMatrixD &aDerivatives) {
	if (theMeasurements.size())
		theMeasurements.back().addGlobals(aLabels, aDerivatives);
}
#endif

/// Define label of point (by GBLTrajectory constructor)
/**
 * \param [in] aLabel Label identifying point
 */
void GblPoint::setLabel(unsigned int aLabel) {
	theLabel = aLabel;
}

/// Retrieve label of point
unsigned int GblPoint::getLabel() const {
	return theLabel;
}

/// Define offset for point (by GBLTrajectory constructor)
/**
 * \param [in] anOffset Offset number
 */
void GblPoint::setOffset(int anOffset) {
	theOffset = anOffset;
}

/// Retrieve offset for point
int GblPoint::getOffset() const {
	return theOffset;
}

/// Define type for point (by GBLTrajectory constructor)
/**
 * \param [in] aType Type (-1: first, 0: inner, 1: last)
 */
void GblPoint::setType(int aType) {
	theType = aType;
}

/// Retrieve point-to-(previous)point jacobian
const Matrix5d& GblPoint::getP2pJacobian() const {
	return p2pJacobian;
}

/// Define jacobian to previous scatterer (by GBLTrajectory constructor)
/**
 * \param [in] aJac Jacobian
 */
void GblPoint::addPrevJacobian(const Matrix5d &aJac) {
	// to optimize: need only two last rows of inverse
	//  prevJacobian = aJac.inverse();
	//  block matrix algebra
	Matrix23d CA = aJac.block<2, 3>(3, 0) * aJac.block<3, 3>(0, 0).inverse(); // C*A^-1
	Matrix2d DCAB = aJac.block<2, 2>(3, 3) - CA * aJac.block<3, 2>(0, 3); // D - C*A^-1 *B
	Matrix2d DCABInv = DCAB.inverse();
	prevJacobian.block<2, 2>(3, 3) = DCABInv;
	prevJacobian.block<2, 3>(3, 0) = -DCABInv * CA;
}

/// Define jacobian to next scatterer (by GBLTrajectory constructor)
/**
 * \param [in] aJac Jacobian
 */
void GblPoint::addNextJacobian(const Matrix5d &aJac) {
	nextJacobian = aJac;
}

/// Retrieve derivatives of local track model
/**
 * Linearized track model: F_u(q/p,u',u) = J*u + S*u' + d*q/p,
 * W is inverse of S, negated for backward propagation.
 * \param [in] aDirection Propagation direction (>0 forward, else backward)
 * \param [out] matW W
 * \param [out] matWJ W*J
 * \param [out] vecWd W*d
 * \exception std::overflow_error : matrix S is singular.
 */
void GblPoint::getDerivatives(int aDirection, Matrix2d &matW, Matrix2d &matWJ,
		Vector2d &vecWd) const {

	Matrix2d matJ;
	Vector2d vecd;
	if (aDirection < 1) {
		matJ = prevJacobian.block<2, 2>(3, 3);
		matW = -prevJacobian.block<2, 2>(3, 1);
		vecd = prevJacobian.block<2, 1>(3, 0);
	} else {
		matJ = nextJacobian.block<2, 2>(3, 3);
		matW = nextJacobian.block<2, 2>(3, 1);
		vecd = nextJacobian.block<2, 1>(3, 0);
	}

	if (!matW.determinant()) {
		std::cout << " GblPoint::getDerivatives failed to invert matrix "
				<< std::endl;
		std::cout
				<< " Possible reason for singular matrix: multiple GblPoints at same arc-length"
				<< std::endl;
		throw std::overflow_error("Singular matrix inversion exception");
	}
	matW = matW.inverse().eval();
	matWJ = matW * matJ;
	vecWd = matW * vecd;
}

/// Print GblPoint
/**
 * \param [in] level print level (0: minimum, >0: more)
 */
void GblPoint::printPoint(unsigned int level) const {
	std::cout << " GblPoint";
	if (theLabel) {
		std::cout << ", label " << theLabel;
		if (theOffset >= 0) {
			std::cout << ", offset " << theOffset;
		}
	}
	if (theMeasurements.size()) {
		std::cout << ", " << theMeasurements.size() << " measurements";
	}
	std::vector<GblMeasurement>::const_iterator itMeas;
	for (itMeas = theMeasurements.begin(); itMeas < theMeasurements.end();
			++itMeas) {
		itMeas->printMeasurement(0);
	}
	if (scatDim == 2) {
		std::cout << ", thin scatterer";
	}
	if (scatDim == 4) {
		std::cout << ", thick scatterer";
	}

	std::cout << std::endl;
	if (level > 0) {
		for (itMeas = theMeasurements.begin(); itMeas < theMeasurements.end();
				++itMeas) {
			itMeas->printMeasurement(level);
		}
		IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
		if (scatDim) {
			std::cout << "  Scatterer" << std::endl;
			std::cout << "   Residuals: "
					<< scatResiduals.head(scatDim).transpose().format(CleanFmt)
					<< std::endl;
			std::cout << "   Precision: "
					<< scatPrecision.head(scatDim).transpose().format(CleanFmt)
					<< std::endl;
		}
		std::cout << "  Jacobian " << std::endl;
		std::cout << "   Point-to-point " << std::endl
				<< p2pJacobian.format(CleanFmt) << std::endl;
		if (theLabel) {
			if (!isFirst())
				std::cout << "   To previous offset (last 2 rows)" << std::endl
						<< prevJacobian.bottomRows<2>().format(CleanFmt)
						<< std::endl;
			if (!isLast())
				std::cout << "   To next offset " << std::endl
						<< nextJacobian.format(CleanFmt) << std::endl;
		}
	}
}

/// Get GblMeasurement iterator for begin
std::vector<GblMeasurement>::iterator GblPoint::getMeasBegin() {
	return theMeasurements.begin();
}

/// Get GblMeasurement iterator for end
std::vector<GblMeasurement>::iterator GblPoint::getMeasEnd() {
	return theMeasurements.end();
}

/// Retrieve global derivatives from a measurement at a point for a single row.
/**
 * \param [in] aMeas  Measurement number
 * \param [in] aRow  Row number
 * \param [out] aLabels Global labels
 * \param [out] aDerivatives  Global derivatives
 */
void GblPoint::getGlobalLabelsAndDerivatives(unsigned int aMeas,
		unsigned int aRow, std::vector<int> &aLabels,
		std::vector<double> &aDerivatives) const {
	theMeasurements[aMeas].getGlobalLabelsAndDerivatives(aRow, aLabels,
			aDerivatives);
}

/// Is first point on (sub)trajectory?
bool GblPoint::isFirst() const {
	return (theType < 0);
}

/// Is last point on (sub)trajectory?
bool GblPoint::isLast() const {
	return (theType > 0);
}

}
