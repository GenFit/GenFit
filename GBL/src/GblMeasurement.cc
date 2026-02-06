/*
 * GBLMeasurement.cpp
 *
 *  Created on: 31 Mar 2021
 *      Author: kleinwrt
 */

/** \file
 *  GblMeasurement methods.
 *
 *  \author Claus Kleinwort, DESY, 2021 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2021 - 2023 Deutsches Elektronen-Synchroton,
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

#include "GblMeasurement.h"
using namespace Eigen;

//! Namespace for the general broken lines package
namespace gbl {

/// Create a measurement.
/**
 * Create measurement (in meas. system) with diagonal or arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aProjection Projection from local to measurement system
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (vector (with diagonal) or (full) matrix)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
GblMeasurement::GblMeasurement(const Eigen::MatrixXd &aProjection,
		const Eigen::VectorXd &aResiduals, const Eigen::MatrixXd &aPrecision,
		double minPrecision) {
	enabled = true;
	measDim = aResiduals.rows();
	measPrecMin = minPrecision;
	if (aPrecision.cols() > 1) {
		// arbitrary precision matrix
		SelfAdjointEigenSolver<MatrixXd> measEigen(aPrecision);
		measTransformation = measEigen.eigenvectors();
		measTransformation.transposeInPlace();
		transFlag = true;
		measResiduals.tail(measDim) = measTransformation * aResiduals;
		measPrecision.tail(measDim) = measEigen.eigenvalues();
		measProjection.bottomRightCorner(measDim, measDim) = measTransformation
				* aProjection;
	} else {
		// diagonal precision matrix
		transFlag = false;
		measResiduals.tail(measDim) = aResiduals;
		measPrecision.tail(measDim) = aPrecision;
		measProjection.bottomRightCorner(measDim, measDim) = aProjection;
	}
}

/// Create a measurement.
/**
 * Create measurement in local system with diagonal or arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (vector (with diagonal) or (full) matrix)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
GblMeasurement::GblMeasurement(const Eigen::VectorXd &aResiduals,
		const Eigen::MatrixXd &aPrecision, double minPrecision) {
	enabled = true;
	measDim = aResiduals.rows();
	measPrecMin = minPrecision;
	if (aPrecision.cols() > 1) {
		// arbitrary precision matrix
		SelfAdjointEigenSolver<MatrixXd> measEigen(aPrecision);
		measTransformation = measEigen.eigenvectors();
		measTransformation.transposeInPlace();
		transFlag = true;
		measResiduals.tail(measDim) = measTransformation * aResiduals;
		measPrecision.tail(measDim) = measEigen.eigenvalues();
		measProjection.bottomRightCorner(measDim, measDim) = measTransformation;
	} else {
		// diagonal precision matrix
		transFlag = false;
		measResiduals.tail(measDim) = aResiduals;
		measPrecision.tail(measDim) = aPrecision;
		measProjection.setIdentity();
	}
}

GblMeasurement::~GblMeasurement() {
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Create a measurement.
/**
 * Create measurement (in meas. system) with diagonal precision (inverse covariance) matrix.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aProjection Projection from local to measurement system
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (diagonal)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
GblMeasurement::GblMeasurement(const TMatrixD &aProjection,
		const TVectorD &aResiduals, const TVectorD &aPrecision,
		double minPrecision) {
	enabled = true;
	measDim = aResiduals.GetNrows();
	measPrecMin = minPrecision;
	transFlag = false;
	unsigned int iOff = 5 - measDim;
	for (unsigned int i = 0; i < measDim; ++i) {
		measResiduals(iOff + i) = aResiduals[i];
		measPrecision(iOff + i) = aPrecision[i];
		for (unsigned int j = 0; j < measDim; ++j) {
			measProjection(iOff + i, iOff + j) = aProjection(i, j);
		}
	}
}

/// Create a measurement.
/**
 * Create measurement (in meas. system) with arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aProjection Projection from local to measurement system
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (matrix)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
GblMeasurement::GblMeasurement(const TMatrixD &aProjection,
		const TVectorD &aResiduals, const TMatrixDSym &aPrecision,
		double minPrecision) {
	enabled = true;
	measDim = aResiduals.GetNrows();
	measPrecMin = minPrecision;
	TMatrixDSymEigen measEigen(aPrecision);
	TMatrixD tmpTransformation(measDim, measDim);
	tmpTransformation = measEigen.GetEigenVectors();
	tmpTransformation.T();
	transFlag = true;
	TVectorD transResiduals = tmpTransformation * aResiduals;
	TVectorD transPrecision = measEigen.GetEigenValues();
	TMatrixD transProjection = tmpTransformation * aProjection;
	measTransformation.resize(measDim, measDim);
	unsigned int iOff = 5 - measDim;
	for (unsigned int i = 0; i < measDim; ++i) {
		measResiduals(iOff + i) = transResiduals[i];
		measPrecision(iOff + i) = transPrecision[i];
		for (unsigned int j = 0; j < measDim; ++j) {
			measTransformation(i, j) = tmpTransformation(i, j);
			measProjection(iOff + i, iOff + j) = transProjection(i, j);
		}
	}
}

/// Create a measurement.
/**
 * Create measurement in local system with diagonal precision (inverse covariance) matrix.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (diagonal)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
GblMeasurement::GblMeasurement(const TVectorD &aResiduals,
		const TVectorD &aPrecision, double minPrecision) {
	enabled = true;
	measDim = aResiduals.GetNrows();
	measPrecMin = minPrecision;
	transFlag = false;
	unsigned int iOff = 5 - measDim;
	for (unsigned int i = 0; i < measDim; ++i) {
		measResiduals(iOff + i) = aResiduals[i];
		measPrecision(iOff + i) = aPrecision[i];
	}
	measProjection.setIdentity();
}

/// Create a measurement.
/**
 * Create measurement in local system with arbitrary precision (inverse covariance) matrix.
 * Will be diagonalized.
 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
 * \param [in] aResiduals Measurement residuals
 * \param [in] aPrecision Measurement precision (matrix)
 * \param [in] minPrecision Minimal precision to accept measurement
 */
GblMeasurement::GblMeasurement(const TVectorD &aResiduals,
		const TMatrixDSym &aPrecision, double minPrecision) {
	enabled = true;
	measDim = aResiduals.GetNrows();
	measPrecMin = minPrecision;
	TMatrixDSymEigen measEigen(aPrecision);
	TMatrixD tmpTransformation(measDim, measDim);
	tmpTransformation = measEigen.GetEigenVectors();
	tmpTransformation.T();
	transFlag = true;
	TVectorD transResiduals = tmpTransformation * aResiduals;
	TVectorD transPrecision = measEigen.GetEigenValues();
	measTransformation.resize(measDim, measDim);
	unsigned int iOff = 5 - measDim;
	for (unsigned int i = 0; i < measDim; ++i) {
		measResiduals(iOff + i) = transResiduals[i];
		measPrecision(iOff + i) = transPrecision[i];
		for (unsigned int j = 0; j < measDim; ++j) {
			measTransformation(i, j) = tmpTransformation(i, j);
			measProjection(iOff + i, iOff + j) = measTransformation(i, j);
		}
	}
}
#endif

/// Add local derivatives.
/**
 * Add local derivatives to measurement.
 * \param [in] aDerivatives Local derivatives (matrix)
 */
void GblMeasurement::addLocals(const Eigen::MatrixXd &aDerivatives) {
	if (measDim) {
		localDerivatives.resize(aDerivatives.rows(), aDerivatives.cols());
		if (transFlag) {
			localDerivatives = measTransformation * aDerivatives;
		} else {
			localDerivatives = aDerivatives;
		}
	}
}

/// Add global derivatives.
/**
 * Add global derivatives to measurement.
 * \param [in] aLabels Global derivatives labels
 * \param [in] aDerivatives Global derivatives (matrix)
 */
void GblMeasurement::addGlobals(const std::vector<int> &aLabels,
		const Eigen::MatrixXd &aDerivatives) {
	if (measDim) {
		globalLabels = aLabels;
		globalDerivatives.resize(aDerivatives.rows(), aDerivatives.cols());
		if (transFlag) {
			globalDerivatives = measTransformation * aDerivatives;
		} else {
			globalDerivatives = aDerivatives;
		}

	}
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Add local derivatives.
/**
 * Add local derivatives to measurement.
 * \param [in] aDerivatives Local derivatives (matrix)
 */
void GblMeasurement::addLocals(const TMatrixD &aDerivatives) {
	if (measDim) {
		unsigned int numDer = aDerivatives.GetNcols();
		localDerivatives.resize(measDim, numDer);
		// convert from ROOT
		MatrixXd tmpDerivatives(measDim, numDer);
		for (unsigned int i = 0; i < measDim; ++i) {
			for (unsigned int j = 0; j < numDer; ++j)
				tmpDerivatives(i, j) = aDerivatives(i, j);
		}
		if (transFlag) {
			localDerivatives = measTransformation * tmpDerivatives;
		} else {
			localDerivatives = tmpDerivatives;
		}
	}
}

/// Add global derivatives.
/**
 * Add global derivatives to measurement.
 * \param [in] aLabels Global derivatives labels
 * \param [in] aDerivatives Global derivatives (matrix)
 */
void GblMeasurement::addGlobals(const std::vector<int> &aLabels,
		const TMatrixD &aDerivatives) {
	if (measDim) {
		globalLabels = aLabels;
		unsigned int numDer = aDerivatives.GetNcols();
		globalDerivatives.resize(measDim, numDer);
		// convert from ROOT
		MatrixXd tmpDerivatives(measDim, numDer);
		for (unsigned int i = 0; i < measDim; ++i) {
			for (unsigned int j = 0; j < numDer; ++j)
				tmpDerivatives(i, j) = aDerivatives(i, j);
		}
		if (transFlag) {
			globalDerivatives = measTransformation * tmpDerivatives;
		} else {
			globalDerivatives = tmpDerivatives;
		}
	}
}
#endif

/// Set enabled flag.
/**
 * Allows to disable measurements at points (in input list of points for trajectory constructor)
 * to construct new trajectory with reduced number of measurements (e.g. after solution of ambiguities).
 *
 * \param [in] flag  enabled flag
 */
void GblMeasurement::setEnabled(bool flag) {
	enabled = flag;
}

/// Get enabled flag.
/**
 * \return enabled flag
 */
bool GblMeasurement::isEnabled() const {
	return enabled;
}

/// Get measurement dimension.
/**
 * Get dimension of measurement (0 = none).
 * \return measurement dimension
 */
unsigned int GblMeasurement::getMeasDim() const {
	return measDim;
}

/// get precision cutoff.
/**
 * \return minimal measurement precision (for usage)
 */
double GblMeasurement::getMeasPrecMin() const {
	return measPrecMin;
}

/// Retrieve measurement of a point.
/**
 * \param [out] aProjection Projection from (diagonalized) measurement to local system
 * \param [out] aResiduals Measurement residuals
 * \param [out] aPrecision Measurement precision (diagonal)
 */
void GblMeasurement::getMeasurement(Matrix5d &aProjection, Vector5d &aResiduals,
		Vector5d &aPrecision) const {
	aProjection.bottomRightCorner(measDim, measDim) =
			measProjection.bottomRightCorner(measDim, measDim);
	aResiduals.tail(measDim) = measResiduals.tail(measDim);
	aPrecision.tail(measDim) = measPrecision.tail(measDim);
}

/// Get measurement transformation (from diagonalization).
/**
 * \param [out] aTransformation Transformation matrix
 */
void GblMeasurement::getMeasTransformation(MatrixXd &aTransformation) const {
	aTransformation.resize(measDim, measDim);
	if (transFlag) {
		aTransformation = measTransformation;
	} else {
		aTransformation.setIdentity();
	}
}

/// Retrieve number of local derivatives from a point.
unsigned int GblMeasurement::getNumLocals() const {
	return localDerivatives.cols();
}

/// Retrieve local derivatives from a point.
const MatrixXd& GblMeasurement::getLocalDerivatives() const {
	return localDerivatives;
}

/// Retrieve number of global derivatives from a point.
unsigned int GblMeasurement::getNumGlobals() const {
	return globalDerivatives.cols();
}

/// Retrieve global derivatives labels from a point.
/**
 * \param [out] aLabels Global labels
 */
void GblMeasurement::getGlobalLabels(std::vector<int> &aLabels) const {
	aLabels = globalLabels;
}

/// Retrieve global derivatives from a point.
/**
 * \param [out] aDerivatives  Global derivatives
 */
void GblMeasurement::getGlobalDerivatives(MatrixXd &aDerivatives) const {
	aDerivatives = globalDerivatives;
}

/// Retrieve global derivatives from a point for a single row.
/**
 * \param [in] aRow  Row number
 * \param [out] aLabels Global labels
 * \param [out] aDerivatives  Global derivatives
 */
void GblMeasurement::getGlobalLabelsAndDerivatives(unsigned int aRow,
		std::vector<int> &aLabels, std::vector<double> &aDerivatives) const {
	aLabels.resize(globalDerivatives.cols());
	aDerivatives.resize(globalDerivatives.cols());
	for (unsigned int i = 0; i < globalDerivatives.cols(); ++i) {
		aLabels[i] = globalLabels[i];
		aDerivatives[i] = globalDerivatives(aRow, i);
	}
}

/// Print GblMeasurement
/**
 * \param [in] level print level (0: minimum, >0: more)
 */
void GblMeasurement::printMeasurement(unsigned int level) const {
	if (level > 0)
		std::cout << "  GblMeasurement";

	if (measDim) {
		std::cout << ", " << measDim << " dimensions";
	}
	if (!enabled) {
		std::cout << ", disabled";
	}
	if (transFlag) {
		std::cout << ", diagonalized";
	}
	if (localDerivatives.cols()) {
		std::cout << ", " << localDerivatives.cols() << " local derivatives";
	}
	if (globalDerivatives.cols()) {
		std::cout << ", " << globalDerivatives.cols() << " global derivatives";
	}
	if (level > 0) {
		std::cout << std::endl;
		IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
		if (measDim) {
			std::cout << "   Projection: " << std::endl
					<< measProjection.bottomRightCorner(measDim, measDim).format(
							CleanFmt) << std::endl;
			std::cout << "   Residuals: "
					<< measResiduals.tail(measDim).transpose().format(CleanFmt)
					<< std::endl;
			std::cout << "   Precision (min.: " << measPrecMin << "): "
					<< measPrecision.tail(measDim).transpose().format(CleanFmt)
					<< std::endl;
		}
		if (localDerivatives.cols()) {
			std::cout << "  Local Derivatives:" << std::endl
					<< localDerivatives.format(CleanFmt) << std::endl;
		}
		if (globalDerivatives.cols()) {
			std::cout << "  Global Labels:";
			for (unsigned int i = 0; i < globalLabels.size(); ++i) {
				std::cout << " " << globalLabels[i];
			}
			std::cout << std::endl;
			std::cout << "  Global Derivatives:"
					<< globalDerivatives.format(CleanFmt) << std::endl;
		}
	}
}

}
