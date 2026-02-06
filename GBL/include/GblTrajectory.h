/*
 * GblTrajectory.h
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \file
 *  GblTrajectory definition.
 *
 *  \author Claus Kleinwort, DESY, 2011 (Claus.Kleinwort@desy.de)
 *  \author Gregor Mittag, DESY, 2017 (templates and other optimizations)
 *
 *
 *  \copyright
 *  Copyright (c) 2011 - 2024 Deutsches Elektronen-Synchroton,
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

#ifndef GBLTRAJECTORY_H_
#define GBLTRAJECTORY_H_

#include <array>
#include "GblPoint.h"
#include "GblData.h"
#include "GblPoint.h"
#include "BorderedBandMatrix.h"
#include "MilleBinary.h"

//! Namespace for the general broken lines package
namespace gbl {

/// GBL trajectory.
/**
 * List of GblPoints ordered by arc length.
 * Can be fitted and optionally written to MP-II binary file.
 */
class GblTrajectory {
public:
	GblTrajectory(const std::vector<GblPoint> &aPointList, bool flagCurv = true,
			bool flagU1dir = true, bool flagU2dir = true);
	GblTrajectory(
			const std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > &aPointsAndTransList);

	/// Create new (simple) trajectory from list of points with external seed.
	/**
	 * Curved trajectory in space (default) or without curvature (q/p) or in one
	 * plane (u-direction) only.
	 * \tparam Seed   Seed precision matrix
	 * \param [in] aPointList List of points
	 * \param [in] aLabel (Signed) label of point for external seed
	 * (<0: in front, >0: after point, slope changes at scatterer!)
	 * \param [in] aSeed Precision matrix of external seed
	 * \param [in] flagCurv Use q/p
	 * \param [in] flagU1dir Use in u1 direction
	 * \param [in] flagU2dir Use in u2 direction
	 */
	template<typename Seed>
	GblTrajectory(const std::vector<GblPoint> &aPointList, unsigned int aLabel,
			const Eigen::MatrixBase<Seed> &aSeed, bool flagCurv = true,
			bool flagU1dir = true, bool flagU2dir = true);

	/// Create new composed trajectory from list of points and transformations with arbitrary external measurements.
	/**
	 * Composed of curved trajectories in space. The precision matrix for the external measurements is specified as matrix.
	 *
	 * \tparam Derivatives  External derivatives
	 * \tparam Measurements   Residuals vector
	 * \tparam Precision   Precision matrix or vector (with diagonal)
	 * \param [in] aPointsAndTransList List containing pairs with list of points and transformation (at inner (first) point)
	 * \param [in] extDerivatives Derivatives of external measurements vs external parameters
	 * \param [in] extMeasurements External measurements (residuals)
	 * \param [in] extPrecisions Precision of external measurements (matrix)
	 */
	template<typename Derivatives, typename Measurements, typename Precision,
			typename std::enable_if<(Precision::ColsAtCompileTime != 1)>::type* =
					nullptr>
	GblTrajectory(
			const std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > &aPointsAndTransList,
			const Eigen::MatrixBase<Derivatives> &extDerivatives,
			const Eigen::MatrixBase<Measurements> &extMeasurements,
			const Eigen::MatrixBase<Precision> &extPrecisions);

	/// Create new composed trajectory from list of points and transformations with independent external measurements.
	/**
	 * Composed of curved trajectories in space. The (diagonal) precision matrix for the external measurements is specified as vector
	 * (containing the diagonal).
	 *
	 * \tparam Derivatives  External derivatives
	 * \tparam Measurements   Residuals vector
	 * \tparam Precision   Precision matrix or vector (with diagonal)
	 * \param [in] aPointsAndTransList List containing pairs with list of points and transformation (at inner (first) point)
	 * \param [in] extDerivatives Derivatives of external measurements vs external parameters
	 * \param [in] extMeasurements External measurements (residuals)
	 * \param [in] extPrecisions Precision of external measurements (vector with diagonal)
	 */
	template<typename Derivatives, typename Measurements, typename Precision,
			typename std::enable_if<(Precision::ColsAtCompileTime == 1)>::type* =
					nullptr>
	GblTrajectory(
			const std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > &aPointsAndTransList,
			const Eigen::MatrixBase<Derivatives> &extDerivatives,
			const Eigen::MatrixBase<Measurements> &extMeasurements,
			const Eigen::MatrixBase<Precision> &extPrecisions);

#ifdef GBL_EIGEN_SUPPORT_ROOT
	// input from ROOT
	GblTrajectory(const std::vector<GblPoint> &aPointList, unsigned int aLabel,
			const TMatrixDSym &aSeed, bool flagCurv = true, bool flagU1dir =
					true, bool flagU2dir = true);
	GblTrajectory(
			const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> > &aPointsAndTransList);
	GblTrajectory(
			const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> > &aPointsAndTransList,
			const TMatrixD &extDerivatives, const TVectorD &extMeasurements,
			const TVectorD &extPrecisions);
	GblTrajectory(
			const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> > &aPointsAndTransList,
			const TMatrixD &extDerivatives, const TVectorD &extMeasurements,
			const TMatrixDSym &extPrecisions);
#endif
	virtual ~GblTrajectory();
	bool isValid() const;
	unsigned int getNumPoints() const;
	unsigned int getExtResults(Eigen::VectorXd &extPar,
			Eigen::MatrixXd &extCov) const;
	unsigned int getResults(int aSignedLabel, Eigen::VectorXd &localPar,
			Eigen::MatrixXd &localCov) const;
	unsigned int getMeasResults(unsigned int aLabel, unsigned int &numData,
			Eigen::VectorXd &aResiduals, Eigen::VectorXd &aMeasErrors,
			Eigen::VectorXd &aResErrors, Eigen::VectorXd &aDownWeights);
	unsigned int getMeasResults(unsigned int aLabel, unsigned int &numData,
			Eigen::VectorXd &aResiduals, Eigen::VectorXd &aMeasErrors);
	unsigned int getScatResults(unsigned int aLabel, unsigned int &numData,
			Eigen::VectorXd &aResiduals, Eigen::VectorXd &aMeasErrors,
			Eigen::VectorXd &aResErrors, Eigen::VectorXd &aDownWeights);
#ifdef GBL_EIGEN_SUPPORT_ROOT
	// input from ROOT
	unsigned int getExtResults(TVectorD &extPar,
			TMatrixDSym &extCov) const;
	unsigned int getResults(int aSignedLabel, TVectorD &localPar,
			TMatrixDSym &localCov) const;
	unsigned int getMeasResults(unsigned int aLabel, unsigned int &numData,
			TVectorD &aResiduals, TVectorD &aMeasErrors, TVectorD &aResErrors,
			TVectorD &aDownWeights);
	unsigned int getScatResults(unsigned int aLabel, unsigned int &numData,
			TVectorD &aResiduals, TVectorD &aMeasErrors, TVectorD &aResErrors,
			TVectorD &aDownWeights);
#endif
	unsigned int getLabels(std::vector<unsigned int> &aLabelList) const;
	unsigned int getLabels(
			std::vector<std::vector<unsigned int> > &aLabelList) const;
	unsigned int fit(double &Chi2, int &Ndf, double &lostWeight,
			const std::string &optionList = "", unsigned int aLabel = 0);
	void milleOut(MilleBinary &aMille);
	void printTrajectory(unsigned int level = 0) const;
	void printPoints(unsigned int level = 0) const;
	void printData() const;
	double getBandCondition() const;
	std::vector<GblData> getData() {return theData;} // not part of original GBL

private:
	unsigned int numAllPoints; ///< Number of all points on trajectory
	std::vector<unsigned int> numPoints; ///< Number of points on (sub)trajectory
	unsigned int numTrajectories; ///< Number of trajectories (in composed trajectory)
	unsigned int numOffsetPoints; ///< Number of points with offsets on trajectory
	unsigned int numOffsets; ///< Number of (1D or 2D) offsets on trajectory
	unsigned int numInnerTransformations; ///< Number of inner transformations to external parameters
	unsigned int numInnerTransOffsets; ///< Number of (points with) offsets affected by inner transformations to external parameters
	unsigned int numCurvature; ///< Number of curvature parameters (0 or 1) or external parameters
	unsigned int numParameters; ///< Number of fit parameters
	unsigned int numLocals; ///< Total number of (additional) local parameters
	unsigned int numMeasurements; ///< Total number of measurements
	unsigned int externalPoint; ///< Label of external point (or 0)
	unsigned int skippedMeasLabel; ///< Label of point with measurements skipped in fit (for unbiased residuals) (or 0)
	unsigned int maxNumGlobals; ///< Max. number of global labels/derivatives per point
	bool constructOK; ///< Trajectory has been successfully constructed (ready for fit/output)
	bool fitOK; ///< Trajectory has been successfully fitted (results are valid)
	std::vector<unsigned int> theDimension; ///< List of active dimensions (0=u1, 1=u2) in fit
	std::vector<std::vector<GblPoint> > thePoints; ///< (list of) List of points on trajectory
	std::vector<GblData> theData; ///< List of data blocks
	std::vector<unsigned int> measDataIndex; ///< mapping points to data blocks from measurements
	std::vector<unsigned int> scatDataIndex; ///< mapping points to data blocks from scatterers
	Eigen::MatrixXd externalSeed; ///< Precision (inverse covariance matrix) of external seed
	// composed trajectory
	std::vector<Eigen::MatrixXd> innerTransformations; ///< Transformations at innermost points of composed trajectory (from common external parameters)
	std::vector<Eigen::MatrixXd> innerTransDer; ///< Derivatives at innermost points of composed trajectory
	std::vector<std::array<unsigned int, 5> > innerTransLab; ///< Labels at innermost points of composed trajectory
	Eigen::MatrixXd externalDerivatives; ///< Derivatives for external measurements of composed trajectory
	Eigen::VectorXd externalMeasurements; ///< Residuals for external measurements of composed trajectory
	Eigen::VectorXd externalPrecisions; ///< Precisions for external measurements of composed trajectory
	// linear equation system
	VVector theVector; ///< Vector of linear equation system
	BorderedBandMatrix theMatrix; ///< (Bordered band) matrix of linear equation system

	std::pair<std::vector<unsigned int>, Eigen::MatrixXd> getJacobian(
			int aSignedLabel) const;
	void getFitToLocalJacobian(std::array<unsigned int, 5> &anIndex,
			Matrix5d &aJacobian, const GblPoint &aPoint, unsigned int measDim,
			unsigned int nJacobian = 1) const;
	unsigned int getFitToKinkJacobian(std::array<unsigned int, 9> &anIndex,
			Matrix49d &aJacobian, const GblPoint &aPoint) const;
	unsigned int getFitToStepJacobian(std::array<unsigned int, 9> &anIndex,
			Matrix49d &aJacobian, const GblPoint &aPoint) const;
	unsigned int getFitToKinkAndStepJacobian(
			std::array<unsigned int, 9> &anIndex, Matrix49d &aJacobian,
			const GblPoint &aPoint) const;
	void construct();
	void defineOffsets();
	void calcJacobians();
	void prepare();
	void buildLinearEquationSystem();
	void predict();
	double downWeight(unsigned int aMethod);
	void getResAndErr(unsigned int aData, bool used, double &aResidual,
			double &aMeasError, double &aResError, double &aDownWeight);
	void getResAndErr(unsigned int aData, double &aResidual,
			double &aMeasError);
};

template<typename Seed>
GblTrajectory::GblTrajectory(const std::vector<GblPoint> &aPointList,
		unsigned int aLabel, const Eigen::MatrixBase<Seed> &aSeed,
		bool flagCurv, bool flagU1dir, bool flagU2dir) :
		numAllPoints(aPointList.size()), numPoints(), numOffsetPoints(0), numOffsets(
				0), numInnerTransformations(0), numInnerTransOffsets(0), numCurvature(
				flagCurv ? 1 : 0), numParameters(0), numLocals(0), numMeasurements(
				0), externalPoint(aLabel), skippedMeasLabel(0), maxNumGlobals(
				0), theDimension(0), thePoints(), theData(), measDataIndex(), scatDataIndex(), externalSeed(
				aSeed), innerTransformations(), externalDerivatives(), externalMeasurements(), externalPrecisions() {

	if (flagU1dir)
		theDimension.push_back(0);
	if (flagU2dir)
		theDimension.push_back(1);
	// simple (single) trajectory
	thePoints.push_back(aPointList);
	numPoints.push_back(numAllPoints);
	construct(); // construct trajectory
}

template<typename Derivatives, typename Measurements, typename Precision,
		typename std::enable_if<(Precision::ColsAtCompileTime != 1)>::type*>
GblTrajectory::GblTrajectory(
		const std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > &aPointsAndTransList,
		const Eigen::MatrixBase<Derivatives> &extDerivatives,
		const Eigen::MatrixBase<Measurements> &extMeasurements,
		const Eigen::MatrixBase<Precision> &extPrecisions) :
		numAllPoints(), numPoints(), numOffsetPoints(0), numOffsets(0), numInnerTransformations(
				aPointsAndTransList.size()), numParameters(0), numLocals(0), numMeasurements(
				0), externalPoint(0), skippedMeasLabel(0), maxNumGlobals(0), theDimension(
				0), thePoints(), theData(), measDataIndex(), scatDataIndex(), externalSeed(), innerTransformations() {

	static_assert(static_cast<int>(Measurements::ColsAtCompileTime) == 1, "GblTrajectory: cols(Measurements) must be 1 (vector)");
	static_assert(static_cast<int>(Measurements::RowsAtCompileTime) == static_cast<int>(Derivatives::RowsAtCompileTime), "GblTrajectory: rows(Measurements) and rows(Derivatives) must be equal");
	static_assert(static_cast<int>(Measurements::RowsAtCompileTime) == static_cast<int>(Precision::RowsAtCompileTime), "GblTrajectory: rows(Measurements) and rows(Precision) must be equal");
	static_assert(static_cast<int>(Precision::RowsAtCompileTime) == static_cast<int>(Precision::ColsAtCompileTime), "GblTrajectory: rows(Precision) and cols(Precision) must be equal");
	// diagonalize external measurement
	Eigen::SelfAdjointEigenSolver<typename Precision::PlainObject> extEigen {
			extPrecisions };
	// @TODO   if (extEigen.info() != Success) abort();
	auto extTransformation = extEigen.eigenvectors().transpose();
	externalDerivatives.resize(extDerivatives.rows(), extDerivatives.cols());
	externalDerivatives = extTransformation * extDerivatives;
	externalMeasurements.resize(extMeasurements.size());
	externalMeasurements = extTransformation * extMeasurements;
	externalPrecisions.resize(extMeasurements.size());
	externalPrecisions = extEigen.eigenvalues();

	for (unsigned int iTraj = 0; iTraj < aPointsAndTransList.size(); ++iTraj) {
		thePoints.push_back(aPointsAndTransList[iTraj].first);
		numPoints.push_back(thePoints.back().size());
		numAllPoints += numPoints.back();
		innerTransformations.push_back(aPointsAndTransList[iTraj].second);
	}
	theDimension.push_back(0);
	theDimension.push_back(1);
	// kinematic (2) or geometric (1) constraint
	numInnerTransOffsets = innerTransformations[0].rows() == 5 ? 2 : 1;
	numCurvature =
			innerTransformations[0].rows() == 5 ?
					innerTransformations[0].cols() :
					innerTransformations[0].cols() + numInnerTransformations;
	construct(); // construct (composed) trajectory
}

template<typename Derivatives, typename Measurements, typename Precision,
		typename std::enable_if<(Precision::ColsAtCompileTime == 1)>::type*>
GblTrajectory::GblTrajectory(
		const std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > &aPointsAndTransList,
		const Eigen::MatrixBase<Derivatives> &extDerivatives,
		const Eigen::MatrixBase<Measurements> &extMeasurements,
		const Eigen::MatrixBase<Precision> &extPrecisions) :
		numAllPoints(), numPoints(), numOffsetPoints(0), numOffsets(0), numInnerTransformations(
				aPointsAndTransList.size()), numParameters(0), numLocals(0), numMeasurements(
				0), externalPoint(0), skippedMeasLabel(0), maxNumGlobals(0), theDimension(
				0), thePoints(), theData(), measDataIndex(), scatDataIndex(), externalSeed(), innerTransformations() {
	static_assert(static_cast<int>(Measurements::ColsAtCompileTime) == 1, "GblTrajectory: cols(Measurements) must be 1 (vector)");
	static_assert(static_cast<int>(Measurements::RowsAtCompileTime) == static_cast<int>(Derivatives::RowsAtCompileTime), "GblTrajectory: rows(Measurements) and rows(Derivatives) must be equal");
	static_assert(static_cast<int>(Measurements::RowsAtCompileTime) == static_cast<int>(Precision::RowsAtCompileTime), "GblTrajectory: rows(Measurements) and rows(Precision) must be equal");

	externalDerivatives = extDerivatives;
	externalMeasurements = extMeasurements;
	externalPrecisions = extPrecisions;

	for (unsigned int iTraj = 0; iTraj < aPointsAndTransList.size(); ++iTraj) {
		thePoints.push_back(aPointsAndTransList[iTraj].first);
		numPoints.push_back(thePoints.back().size());
		numAllPoints += numPoints.back();
		innerTransformations.push_back(aPointsAndTransList[iTraj].second);
	}
	theDimension.push_back(0);
	theDimension.push_back(1);
	// kinematic (2) or geometric (1) constraint
	numInnerTransOffsets = innerTransformations[0].rows() == 5 ? 2 : 1;
	numCurvature =
			innerTransformations[0].rows() == 5 ?
					innerTransformations[0].cols() :
					innerTransformations[0].cols() + numInnerTransformations;
	construct(); // construct (composed) trajectory
}

}
#endif /* GBLTRAJECTORY_H_ */
