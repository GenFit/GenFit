/*
 * GblMeasurement.h
 *
 *  Created on: 31 Mar 2021
 *      Author: kleinwrt
 */

/** \file
 *  GblMeasurement (with optional derivatives) definition.
 *
 *  \author Claus Kleinwort, DESY, 2021 (Claus.Kleinwort@desy.de)
 *
 *
 *  \copyright
 *  Copyright (c) 2021 Deutsches Elektronen-Synchroton,
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

#ifndef GBLMEASUREMENT_H_
#define GBLMEASUREMENT_H_

#include<iostream>
#include<vector>
#include <array>
#include<math.h>
#include "VMatrix.h"

#ifdef GBL_EIGEN_SUPPORT_ROOT
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#endif

#include "Eigen/Dense"

namespace gbl {

typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 2, 3> Matrix23d;
typedef Eigen::Matrix<double, 2, 5> Matrix25d;
typedef Eigen::Matrix<double, 3, 2> Matrix32d;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;

/// Measurement at point.
/**
 * User supplied measurement at point on (initial) trajectory.
 *
 * Must have measurement (1D - 5D). May have:
 *
 *   -# Additional local parameters (with derivatives). Fitted together with track parameters.
 *   -# Additional global parameters (with labels and derivatives). Not fitted, only passed
 *      on to (binary) file for fitting with Millepede-II.
 */
class GblMeasurement {
public:
	GblMeasurement(const Eigen::MatrixXd &aProjection,
			const Eigen::VectorXd &aResiduals,
			const Eigen::MatrixXd &aPrecision, double minPrecision = 0.);
	GblMeasurement(const Eigen::VectorXd &aResiduals,
			const Eigen::MatrixXd &aPrecision, double minPrecision = 0.);
	GblMeasurement(const GblMeasurement&) = default;
	GblMeasurement& operator=(const GblMeasurement&) = default;
	GblMeasurement(GblMeasurement&&) = default;
	GblMeasurement& operator=(GblMeasurement&&) = default;
	virtual ~GblMeasurement();
#ifdef GBL_EIGEN_SUPPORT_ROOT
	// input via ROOT
	GblMeasurement(const TMatrixD &aProjection, const TVectorD &aResiduals,
			const TVectorD &aPrecision, double minPrecision = 0.);
	GblMeasurement(const TMatrixD &aProjection, const TVectorD &aResiduals,
			const TMatrixDSym &aPrecision, double minPrecision = 0.);
	GblMeasurement(const TVectorD &aResiduals, const TVectorD &aPrecision,
			double minPrecision = 0.);
	GblMeasurement(const TVectorD &aResiduals, const TMatrixDSym &aPrecision,
			double minPrecision = 0.);
	void addLocals(const TMatrixD &aDerivatives);
	void addGlobals(const std::vector<int> &aLabels,
			const TMatrixD &aDerivatives);
#endif

	// input via Eigen
	void addLocals(const Eigen::MatrixXd &aDerivatives);
	void addGlobals(const std::vector<int> &aLabels,
			const Eigen::MatrixXd &aDerivatives);
	void setEnabled(bool flag);
	bool isEnabled() const;
	unsigned int getMeasDim() const;
	double getMeasPrecMin() const;
	void getMeasurement(Matrix5d &aProjection, Vector5d &aResiduals,
			Vector5d &aPrecision) const;
	void getMeasTransformation(Eigen::MatrixXd &aTransformation) const;
	unsigned int getNumLocals() const;
	const Eigen::MatrixXd& getLocalDerivatives() const;
	unsigned int getNumGlobals() const;
	void printMeasurement(unsigned int level = 0) const;
	void getGlobalLabels(std::vector<int> &aLabels) const;
	void getGlobalDerivatives(Eigen::MatrixXd &aDerivatives) const;
	void getGlobalLabelsAndDerivatives(unsigned int aRow,
			std::vector<int> &aLabels, std::vector<double> &aDerivatives) const;

private:
	bool enabled; ///< Enabled flag (to be used)
	unsigned int measDim; ///< Dimension of measurement (1-5), 0 indicates absence of measurement
	double measPrecMin; ///< Minimal measurement precision (for usage)
	Matrix5d measProjection; ///< Projection from measurement to local system

	Vector5d measResiduals; ///< Measurement residuals
	Vector5d measPrecision; ///< Measurement precision (diagonal of inverse covariance matrix)
	bool transFlag; ///< Transformation exists?
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
			Eigen::ColMajor /* default */, 5, 5> measTransformation; ///< Transformation of diagonalization (of meas. precision matrix)
	Eigen::MatrixXd localDerivatives; ///< Derivatives of measurement vs additional local (fit) parameters
	std::vector<int> globalLabels; ///< Labels of global (MP-II) derivatives
	Eigen::MatrixXd globalDerivatives; ///< Derivatives of measurement vs additional global (MP-II) parameters
};
}

#endif /* GBLMEASUREMENT_H_ */
