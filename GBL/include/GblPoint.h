/*
 * GblPoint.h
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \file
 *  GblPoint definition.
 *
 *  \author Claus Kleinwort, DESY, 2011 (Claus.Kleinwort@desy.de)
 *  \author Gregor Mittag, DESY, 2017 (templates and other optimizations)
 *
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

#ifndef GBLPOINT_H_
#define GBLPOINT_H_

#include "GblMeasurement.h"

#include<iostream>
#include<vector>
#include<math.h>
#include <stdexcept>

#include "Eigen/Dense"

namespace gbl {

typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix<double, 2, 3> Matrix23d;
typedef Eigen::Matrix<double, 2, 5> Matrix25d;
typedef Eigen::Matrix<double, 3, 2> Matrix32d;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 4, 9> Matrix49d;

/// Point on trajectory
/**
 * User supplied point on (initial) trajectory.
 *
 * Must have jacobian for propagation from previous point. May have:
 *
 *   -# Measurement(s) (1D - 5D)
 *   -# Scatterer (thin, 2D kinks)
 */
class GblPoint {
public:
	GblPoint(const Matrix5d &aJacobian, unsigned int numMeasReserve = 0);
	GblPoint(const GblPoint&) = default;
	GblPoint& operator=(const GblPoint&) = default;
	GblPoint(GblPoint&&) = default;
	GblPoint& operator=(GblPoint&&) = default;
	virtual ~GblPoint();
#ifdef GBL_EIGEN_SUPPORT_ROOT
	// input via ROOT
	GblPoint(const TMatrixD &aJacobian);
	void addMeasurement(const TMatrixD &aProjection, const TVectorD &aResiduals,
			const TVectorD &aPrecision, double minPrecision = 0.);
	void addMeasurement(const TMatrixD &aProjection, const TVectorD &aResiduals,
			const TMatrixDSym &aPrecision, double minPrecision = 0.);
	void addMeasurement(const TVectorD &aResiduals, const TVectorD &aPrecision,
			double minPrecision = 0.);
	void addMeasurement(const TVectorD &aResiduals,
			const TMatrixDSym &aPrecision, double minPrecision = 0.);
	void addScatterer(const TVectorD &aResiduals, const TVectorD &aPrecision);
	void addScatterer(const TVectorD &aResiduals,
			const TMatrixDSym &aPrecision);
	void addThickScatterer(const TVectorD &aResiduals,
			const TMatrixDSym &aPrecision);
	void addLocals(const TMatrixD &aDerivatives);
	void addGlobals(const std::vector<int> &aLabels,
			const TMatrixD &aDerivatives);
#endif
	// input via Eigen

	/// Add a measurement to a point.
	/**
	 * Add measurement (in measurement system) with arbitrary precision (inverse covariance) matrix.
	 * Will be diagonalized.
	 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
	 * \tparam Projection  Projection matrix
	 * \tparam Residuals   Residuals vector
	 * \tparam Precision   Precision matrix or vector (with diagonal)
	 * \param [in] aProjection Projection from local to measurement system (derivative of measurement vs local parameters)
	 * \param [in] aResiduals Measurement residuals
	 * \param [in] aPrecision Measurement precision (matrix)
	 * \param [in] minPrecision Minimal precision to accept measurement
	 */
	template<typename Projection, typename Residuals, typename Precision,
			typename std::enable_if<(Precision::ColsAtCompileTime != 1)>::type* =
					nullptr>
	void addMeasurement(const Eigen::MatrixBase<Projection> &aProjection,
			const Eigen::MatrixBase<Residuals> &aResiduals,
			const Eigen::MatrixBase<Precision> &aPrecision,
			double minPrecision = 0.);

	/// Add a measurement to a point.
	/**
	 * Add measurement (in measurement system) with diagonal precision (inverse covariance) matrix.
	 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
	 * \tparam Projection  Projection matrix
	 * \tparam Residuals   Residuals vector
	 * \tparam Precision   Precision matrix or vector (with diagonal)
	 * \param [in] aProjection Projection from local to measurement system (derivative of measurement vs local parameters)
	 * \param [in] aResiduals Measurement residuals
	 * \param [in] aPrecision Measurement precision (vector with diagonal)
	 * \param [in] minPrecision Minimal precision to accept measurement
	 */
	template<typename Projection, typename Residuals, typename Precision,
			typename std::enable_if<(Precision::ColsAtCompileTime == 1)>::type* =
					nullptr>
	void addMeasurement(const Eigen::MatrixBase<Projection> &aProjection,
			const Eigen::MatrixBase<Residuals> &aResiduals,
			const Eigen::MatrixBase<Precision> &aPrecision,
			double minPrecision = 0.);

	/// Add a measurement to a point.
	/**
	 * Add measurement in local system with arbitrary precision (inverse covariance) matrix.
	 * Will be diagonalized.
	 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
	 * \tparam Residuals   Residuals vector
	 * \tparam Precision   Precision matrix or vector (with diagonal)
	 * \param [in] aResiduals Measurement residuals
	 * \param [in] aPrecision Measurement precision (matrix)
	 * \param [in] minPrecision Minimal precision to accept measurement
	 */
	template<typename Residuals, typename Precision, typename std::enable_if<
			(Precision::ColsAtCompileTime != 1)>::type* = nullptr>
	void addMeasurement(const Eigen::MatrixBase<Residuals> &aResiduals,
			const Eigen::MatrixBase<Precision> &aPrecision,
			double minPrecision = 0.);

	/// Add a measurement to a point.
	/**
	 * Add measurement in local system with diagonal precision (inverse covariance) matrix.
	 * ((up to) 2D: position, 4D: slope+position, 5D: curvature+slope+position)
	 * \tparam Residuals   Residuals vector
	 * \tparam Precision   Precision matrix or vector (with diagonal)
	 * \param [in] aResiduals Measurement residuals
	 * \param [in] aPrecision Measurement precision (vector with diagonal)
	 * \param [in] minPrecision Minimal precision to accept measurement
	 */
	template<typename Residuals, typename Precision, typename std::enable_if<
			(Precision::ColsAtCompileTime == 1)>::type* = nullptr>
	void addMeasurement(const Eigen::MatrixBase<Residuals> &aResiduals,
			const Eigen::MatrixBase<Precision> &aPrecision,
			double minPrecision = 0.);

	void addScatterer(const Eigen::Vector2d &aResiduals,
			const Eigen::Matrix2d &aPrecision);
	void addScatterer(const Eigen::Vector2d &aResiduals,
			const Eigen::Vector2d &aPrecision);
	void addThickScatterer(const Eigen::Vector4d &aResiduals,
			const Eigen::Matrix4d &aPrecision);

	/// Add local derivatives to a point.
	/**
	 * Point needs to have a measurement.
	 * \tparam Derivative  Derivatives matrix
	 * \param [in] aDerivatives Local derivatives (matrix)
	 */
	template<typename Derivative>
	void addLocals(const Eigen::MatrixBase<Derivative> &aDerivatives);
	template<typename Derivative>

	/// Add global derivatives to a point.
	/**
	 * Point needs to have a measurement.
	 * \tparam Derivative  Derivatives matrix
	 * \param [in] aLabels Global derivatives labels
	 * \param [in] aDerivatives Global derivatives (matrix)
	 */
	void addGlobals(const std::vector<int> &aLabels,
			const Eigen::MatrixBase<Derivative> &aDerivatives);
	//
	unsigned int numMeasurements() const;
	unsigned int getScatDim() const;
	void getScatterer(Eigen::Matrix4d &aTransformation,
			Eigen::Vector4d &aResiduals, Eigen::Vector4d &aPrecision) const;
	void getReducedScatterer(Eigen::Matrix4d &aTransformation,
			Eigen::Vector4d &aResiduals, Eigen::Vector4d &aPrecision) const;
	void getScatTransformation(Eigen::MatrixXd &aTransformation) const;
	unsigned int getLabel() const;
	int getOffset() const;
	const Matrix5d& getP2pJacobian() const;
	void getDerivatives(int aDirection, Eigen::Matrix2d &matW,
			Eigen::Matrix2d &matWJ, Eigen::Vector2d &vecWd) const;
	void printPoint(unsigned int level = 0) const;
	std::vector<GblMeasurement>::iterator getMeasBegin();
	std::vector<GblMeasurement>::iterator getMeasEnd();
	void getGlobalLabelsAndDerivatives(unsigned int aMeas, unsigned int aRow,
			std::vector<int> &aLabels, std::vector<double> &aDerivatives) const;
	bool isFirst() const;
	bool isLast() const;

private:
	friend class GblTrajectory; // to have the following setters private
	void setLabel(unsigned int aLabel);
	void setOffset(int anOffset);
	void setType(int aType);
	void addPrevJacobian(const Matrix5d &aJac);
	void addNextJacobian(const Matrix5d &aJac);

	unsigned int theLabel; ///< Label identifying point
	int theOffset; ///< Offset number at point if not negative (else interpolation needed)
	int theType; ///< Type (-1: first, 0: inner, 1: last)
	Matrix5d p2pJacobian  = Matrix5d::Zero(5,5); ///< Point-to-point jacobian from previous point
	Matrix5d prevJacobian = Matrix5d::Zero(5,5); ///< Jacobian to previous scatterer (or first measurement)
	Matrix5d nextJacobian = Matrix5d::Zero(5,5); ///< Jacobian to next scatterer (or last measurement)
	unsigned int scatDim; ///< Dimension of scatterer (0: none, 2: thin, 4:thick)
	//Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
	//		Eigen::ColMajor /* default */, 4, 4> scatTransformation; ///< Transformation of diagonalization (of scat. precision matrix)
	Eigen::Matrix4d scatTransformation; ///< Transformation of diagonalization (of scat. precision matrix)
	Eigen::Vector4d scatResiduals; ///< Scattering residuals (initial kinks if iterating)
	Eigen::Vector4d scatPrecision; ///< Scattering precision (diagonal of inverse covariance matrix)
	std::vector<GblMeasurement> theMeasurements; ///< List of measurements at point
};

template<typename Projection, typename Residuals, typename Precision,
		typename std::enable_if<(Precision::ColsAtCompileTime != 1)>::type*>
void GblPoint::addMeasurement(const Eigen::MatrixBase<Projection> &aProjection,
		const Eigen::MatrixBase<Residuals> &aResiduals,
		const Eigen::MatrixBase<Precision> &aPrecision, double minPrecision) {
	static_assert(static_cast<int>(Residuals::ColsAtCompileTime) == 1, "addMeasurement: cols(Residuals) must be 1 (vector)");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) <= 5 or static_cast<int>(Residuals::RowsAtCompileTime) == Eigen::Dynamic, "addMeasurement: rows(Residuals) must be 1-5 or dynamic");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) == static_cast<int>(Precision::RowsAtCompileTime), "addMeasurement: rows(Residuals) and rows(Precision) must be equal");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) == static_cast<int>(Projection::RowsAtCompileTime), "addMeasurement: rows(Residuals) and rows(Projection) must be equal");
	static_assert(static_cast<int>(Precision::RowsAtCompileTime) == static_cast<int>(Precision::ColsAtCompileTime), "addMeasurement: rows(Precision) and cols(Precision) must be equal");
	static_assert(static_cast<int>(Projection::RowsAtCompileTime) == static_cast<int>(Projection::ColsAtCompileTime), "addMeasurement: rows(Projection) and cols(Projection) must be equal");
	theMeasurements.emplace_back(aProjection, aResiduals, aPrecision,
			minPrecision);
}

template<typename Projection, typename Residuals, typename Precision,
		typename std::enable_if<(Precision::ColsAtCompileTime == 1)>::type*>
void GblPoint::addMeasurement(const Eigen::MatrixBase<Projection> &aProjection,
		const Eigen::MatrixBase<Residuals> &aResiduals,
		const Eigen::MatrixBase<Precision> &aPrecision, double minPrecision) {
	static_assert(static_cast<int>(Residuals::ColsAtCompileTime) == 1, "addMeasurement: cols(Residuals) must be 1 (vector)");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) <= 5 or static_cast<int>(Residuals::RowsAtCompileTime) == Eigen::Dynamic, "addMeasurement: rows(Residuals) must be 1-5 or dynamic");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) == static_cast<int>(Precision::RowsAtCompileTime), "addMeasurement: rows(Residuals) and rows(Precision) must be equal");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) == static_cast<int>(Projection::RowsAtCompileTime), "addMeasurement: rows(Residuals) and rows(Projection) must be equal");
	static_assert(static_cast<int>(Projection::RowsAtCompileTime) == static_cast<int>(Projection::ColsAtCompileTime), "addMeasurement: rows(Projection) and cols(Projection) must be equal");
	//theMeasurements.push_back(GblMeasurement(aProjection, aResiduals, aPrecision, minPrecision));
	theMeasurements.emplace_back(aProjection, aResiduals, aPrecision,
			minPrecision);
}

template<typename Residuals, typename Precision, typename std::enable_if<
		(Precision::ColsAtCompileTime != 1)>::type*>
void GblPoint::addMeasurement(const Eigen::MatrixBase<Residuals> &aResiduals,
		const Eigen::MatrixBase<Precision> &aPrecision, double minPrecision) {
	static_assert(static_cast<int>(Residuals::ColsAtCompileTime) == 1, "addMeasurement: cols(Residuals) must be 1 (vector)");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) <= 5 or static_cast<int>(Residuals::RowsAtCompileTime) == Eigen::Dynamic, "addMeasurement: rows(Residuals) must be 1-5 or dynamic");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) == static_cast<int>(Precision::RowsAtCompileTime), "addMeasurement: rows(Residuals) and rows(Precision) must be equal");
	theMeasurements.emplace_back(aResiduals, aPrecision, minPrecision);
}

template<typename Residuals, typename Precision, typename std::enable_if<
		(Precision::ColsAtCompileTime == 1)>::type*>
void GblPoint::addMeasurement(const Eigen::MatrixBase<Residuals> &aResiduals,
		const Eigen::MatrixBase<Precision> &aPrecision, double minPrecision) {
	static_assert(static_cast<int>(Residuals::ColsAtCompileTime) == 1, "addMeasurement: cols(Residuals) must be 1 (vector)");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) <= 5 or static_cast<int>(Residuals::RowsAtCompileTime) == Eigen::Dynamic, "addMeasurement: rows(Residuals) must be 1-5 or dynamic");
	static_assert(static_cast<int>(Residuals::RowsAtCompileTime) == static_cast<int>(Precision::RowsAtCompileTime), "addMeasurement: rows(Residuals) and rows(Precision) must be equal");
	theMeasurements.emplace_back(aResiduals, aPrecision, minPrecision);
}

template<typename Derivative>
void GblPoint::addLocals(const Eigen::MatrixBase<Derivative> &aDerivatives) {
	if (theMeasurements.size())
		theMeasurements.back().addLocals(aDerivatives);
}

template<typename Derivative>
void GblPoint::addGlobals(const std::vector<int> &aLabels,
		const Eigen::MatrixBase<Derivative> &aDerivatives) {
	if (theMeasurements.size())
		theMeasurements.back().addGlobals(aLabels, aDerivatives);
}

}
#endif /* GBLPOINT_H_ */
