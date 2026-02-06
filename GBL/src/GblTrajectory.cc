/*
 * GblTrajectory.cpp
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \file
 *  GblTrajectory methods.
 *
 *  \author Claus Kleinwort, DESY, 2011 (Claus.Kleinwort@desy.de)
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

/** \mainpage General information
 *
 *  \section intro_sec Introduction
 *
 *  For a track with an initial trajectory from a prefit of the
 *  measurements (internal seed) or an external prediction
 *  (external seed) the description of multiple scattering is
 *  added by offsets in a local system. Along the initial
 *  trajectory points are defined with can describe a measurement
 *  or a scatterer or both. Measurements are arbitrary (linear)
 *  functions of the local track parameters at a point (e.g. 2D:
 *  position, 4D: direction+position). Multiple measurements can be
 *  added to a point (and later disabled) to implement
 *  (and resolve) ambiguities.
 *  The refit provides corrections
 *  to the local track parameters (in the local system) and the
 *  corresponding covariance matrix at any of those points.
 *  Non-diagonal covariance matrices of
 *  measurements will be diagonalized internally.
 *  Outliers can be down-weighted by use of M-estimators.
 *  At one point the measurements can be omitted from the refit
 *  to calculate unbiased residuals.
 *
 *  A position measurement is in a plane defined by two directions.
 *  Along one direction the measurement precision may be zero,
 *  defining a 1D measurement in the other direction.
 *
 *  The broken lines trajectory is defined by (2D) offsets at the
 *  first and last point and all points with a scatterer. The
 *  prediction for a measurement is obtained by interpolation of
 *  the enclosing offsets and for triplets of adjacent offsets
 *  kink angles are determined (*thin* scatterer, single offset)
 *  or for quadruplets of adjacent offsets kink angles and steps
 *  are determined (*thick* scatterer, two offsets, see logo).
 *  This requires for all points the jacobians for propagation
 *  to the previous and next point with an offset.
 *  These are calculated from the point-to-point jacobians along
 *  the initial trajectory. The sequence of points has to be
 *  strictly monotonic in arc-length.
 *
 *  Additional local or global parameters can be added and the
 *  trajectories can be written to special binary files for
 *  calibration and alignment with Millepede-II.
 *  (V. Blobel, NIM A, 566 (2006), pp. 5-13).
 *
 *  Besides simple trajectories describing the path of a single
 *  particle composed trajectories are supported. These are
 *  constructed from the trajectories of multiple particles and
 *  some external parameters (like those describing a decay)
 *  and transformations at the first points from the external
 *  to the local track parameters.
 *  The external parameters can describe the full kinematics (all track parameters)
 *  or only the geometry (position parameters at common vertex).
 *
 *  The conventions for the coordinate systems follow:
 *  Derivation of Jacobians for the propagation of covariance
 *  matrices of track parameters in homogeneous magnetic fields
 *  A. Strandlie, W. Wittek, NIM A, 566 (2006) 687-698.
 *
 *  The source code is available at the DESY GitLab server, see:
 *  https://gitlab.desy.de/claus.kleinwort/general-broken-lines/-/wikis/home
 *
 *  \section call_sec Calling sequence
 *
 *    -# Create list of points on initial trajectory:\n
 *            <tt>std::vector<GblPoint> list</tt>
 *    -# For all points on initial trajectory:
 *        - Create points and add appropriate attributes:\n
 *           - <tt>point = gbl::GblPoint(..)</tt>
 *           - <tt>point.addMeasurement(..)</tt>
 *           - Add additional local or global parameters to measurement:\n
 *             - <tt>point.addLocals(..)</tt>
 *             - <tt>point.addGlobals(..)</tt>
 *           - <tt>point.addScatterer(..)</tt> or <tt>point.addThickScatterer(..)</tt>
 *        - Add point (ordered by arc length) to list:\n
 *            <tt>list.push_back(point)</tt>
 *    -# Create (simple) trajectory from list of points:\n
 *            <tt>traj = gbl::GblTrajectory (list)</tt>
 *    -# Optionally with external seed:\n
 *            <tt>traj = gbl::GblTrajectory (list,seed)</tt>
 *    -# Optionally check validity of trajectory:\n
 *            <tt>if (!traj.isValid()) .. //abort</tt>
 *    -# Fit trajectory (potentially several times with different options), return error code,
 *       get Chi2, Ndf (and weight lost by M-estimators):\n
 *            <tt>ierr = traj.fit(..)</tt>
 *    -# For any point on initial trajectory:
 *        - Get corrections and covariance matrix for track parameters:\n
 *            <tt>[..] = traj.getResults(label)</tt>
 *        - Optionally get residuals with errors for measurements:\n
 *            <tt>[..] = traj.getMeasResults(label) </tt>
 *        - Optionally get residuals with errors for scatterers:\n
 *            <tt>[..] = traj.getScatResults(label) </tt>
 *    -# Optionally write trajectory to MP binary file (doesn't needs to be fitted):\n
 *            <tt>traj.milleOut(..)</tt>
 *
 *  \section loc_sec Local system and local parameters
 *  At each point on the trajectory a local coordinate system with local track
 *  parameters has to be defined. The first of the five parameters describes
 *  the bending, the next two the direction and the last two the position (offsets).
 *  The curvilinear system (T,U,V) with parameters (q/p, lambda, phi, x_t, y_t)
 *  is well suited.
 *
 *  \section impl_sec Implementation
 *
 *  Matrices are implemented with the EIGEN template library (eigen.tuxfamily.org). User input or output is in the
 *  form of MatrixXd. With the preprocessor directive <tt>GBL_EIGEN_SUPPORT_ROOT</tt> user input and output with
 *  ROOT (root.cern.ch) matrices is supported too.
 *  Internally Matrix<n>d are used for fixed sized and simple matrices
 *  based on std::vector<> for variable sized matrices.
 *  Several GBL methods are implemented as templates to allow more EIGEN compile time optimization.
 *
 *  \section jna_sec Java Native Access wrappers
 *  JavaNativeAccessWrappers.cpp are provided by Tom Eichlersmith (Univ. Minnesota).
 *
 *  \section example_sec Examples
 *  Technical examples are given in example1.cpp, example2.cpp, example3.cpp and example4.cpp.
 *  An example silicon tracker is described in exampleSit.cpp
 *  and an example sector of forward drift chambers in exampleDc.cpp.
 *  Composed trajectories in a cylindrical drift chamber are used in exampleComposedGeo.cpp and exampleComposedKin.cpp.
 *
 *  \section ref_sec References
 *    - V. Blobel, C. Kleinwort, F. Meier,
 *      Fast alignment of a complex tracking detector using advanced track models,
 *      Computer Phys. Communications (2011), doi:10.1016/j.cpc.2011.03.017
 *    - C. Kleinwort, General Broken Lines as advanced track fitting method,
 *      NIM A, 673 (2012), 107-110, doi:10.1016/j.nima.2012.01.024
 */

#include "GblTrajectory.h"
using namespace Eigen;

//! Namespace for the general broken lines package
namespace gbl {

/// Create new (simple) trajectory from list of points.
/**
 * Curved trajectory in space (default) or without curvature (q/p) or in one
 * plane (u-direction) only.
 * \param [in] aPointList List of points
 * \param [in] flagCurv Use q/p
 * \param [in] flagU1dir Use in u1 direction
 * \param [in] flagU2dir Use in u2 direction
 */
GblTrajectory::GblTrajectory(const std::vector<GblPoint> &aPointList,
		bool flagCurv, bool flagU1dir, bool flagU2dir) :
		numAllPoints(aPointList.size()), numPoints(), numOffsetPoints(0), numOffsets(0),
				numInnerTransformations(0), numInnerTransOffsets(0), numCurvature(flagCurv ? 1 : 0),
				numParameters(0), numLocals(0), numMeasurements(0),
				externalPoint(0), skippedMeasLabel(0), maxNumGlobals(0), theDimension(0),
				thePoints(), theData(), measDataIndex(), scatDataIndex(), externalSeed(),
				innerTransformations(), externalDerivatives(), externalMeasurements(), externalPrecisions() {

	if (flagU1dir)
		theDimension.push_back(0);
	if (flagU2dir)
		theDimension.push_back(1);
	// simple (single) trajectory
	thePoints.emplace_back(std::move(aPointList));
	numPoints.push_back(numAllPoints);
	construct(); // construct trajectory
}

/// Create new composed trajectory from list of points and transformations.
/**
 * Composed of curved trajectories in space.
 * \param [in] aPointsAndTransList List containing pairs with list of points and transformation (at inner (first) point)
 */
GblTrajectory::GblTrajectory(
		const std::vector<std::pair<std::vector<GblPoint>, Eigen::MatrixXd> > &aPointsAndTransList) :
		numAllPoints(), numPoints(), numOffsetPoints(0), numOffsets(0), numInnerTransformations(
				aPointsAndTransList.size()), numParameters(0), numLocals(0), numMeasurements(
				0), externalPoint(0), skippedMeasLabel(0), maxNumGlobals(0), theDimension(
				0), thePoints(), theData(), measDataIndex(), scatDataIndex(), externalSeed(), innerTransformations(), externalDerivatives(), externalMeasurements(), externalPrecisions() {

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

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Create new (simple) trajectory from list of points with external seed.
/**
 * Curved trajectory in space (default) or without curvature (q/p) or in one
 * plane (u-direction) only.
 * \param [in] aPointList List of points
 * \param [in] aLabel (Signed) label of point for external seed
 * (<0: in front, >0: after point, slope changes at scatterer!)
 * \param [in] aSeed Precision matrix of external seed
 * \param [in] flagCurv Use q/p
 * \param [in] flagU1dir Use in u1 direction
 * \param [in] flagU2dir Use in u2 direction
 */
GblTrajectory::GblTrajectory(const std::vector<GblPoint> &aPointList,
		unsigned int aLabel, const TMatrixDSym &aSeed, bool flagCurv,
		bool flagU1dir, bool flagU2dir) :
numAllPoints(aPointList.size()), numPoints(), numOffsetPoints(0), numOffsets(0),
numInnerTransformations(0), numInnerTransOffsets(0), numCurvature(flagCurv ? 1 : 0), numParameters(0),
numLocals(0), numMeasurements(0), externalPoint(aLabel), skippedMeasLabel(0),
maxNumGlobals(0), theDimension(0), thePoints(), theData(), measDataIndex(),
scatDataIndex(), externalSeed(), innerTransformations(), externalDerivatives(),
externalMeasurements(), externalPrecisions() {

	if (flagU1dir)
	theDimension.push_back(0);
	if (flagU2dir)
	theDimension.push_back(1);
	// convert from ROOT
	unsigned int nParSeed = aSeed.GetNrows();
	externalSeed.resize(nParSeed, nParSeed);
	for (unsigned int i = 0; i < nParSeed; ++i) {
		for (unsigned int j = 0; j < nParSeed; ++j) {
			externalSeed(i, j) = aSeed(i, j);
		}
	}
	// simple (single) trajectory
	thePoints.emplace_back(std::move(aPointList));
	numPoints.push_back(numAllPoints);
	construct();// construct trajectory
}

/// Create new composed trajectory from list of points and transformations.
/**
 * Composed of curved trajectories in space.
 * \param [in] aPointsAndTransList List containing pairs with list of points and transformation (at inner (first) point)
 */
GblTrajectory::GblTrajectory(const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> > &aPointsAndTransList) :
numAllPoints(), numPoints(), numOffsetPoints(0), numOffsets(0), numInnerTransformations(aPointsAndTransList.size()),
numParameters(0), numLocals(0), numMeasurements(0), externalPoint(0),
skippedMeasLabel(0), maxNumGlobals(0), theDimension(0), thePoints(), theData(),
measDataIndex(), scatDataIndex(), externalSeed(), innerTransformations(),
externalDerivatives(), externalMeasurements(), externalPrecisions() {

	for (unsigned int iTraj = 0; iTraj < aPointsAndTransList.size(); ++iTraj) {
		thePoints.emplace_back(std::move(aPointsAndTransList[iTraj].first));
		numPoints.push_back(thePoints.back().size());
		numAllPoints += numPoints.back();
		// convert from ROOT
		unsigned int nRowTrans = aPointsAndTransList[iTraj].second.GetNrows();
		unsigned int nColTrans = aPointsAndTransList[iTraj].second.GetNcols();
		MatrixXd aTrans(nRowTrans, nColTrans);
		for (unsigned int i = 0; i < nRowTrans; ++i) {
			for (unsigned int j = 0; j < nColTrans; ++j) {
				aTrans(i, j) = aPointsAndTransList[iTraj].second(i, j);
			}
		}
		innerTransformations.emplace_back(std::move(aTrans));
	}
	theDimension.push_back(0);
	theDimension.push_back(1);
	// kinematic (2) or geometric (1) constraint
	numInnerTransOffsets = innerTransformations[0].rows() == 5 ? 2 : 1;
	numCurvature =
	innerTransformations[0].rows() == 5 ?
	innerTransformations[0].cols() :
	innerTransformations[0].cols() + numInnerTransformations;
	construct();// construct (composed) trajectory
}

/// Create new composed trajectory from list of points and transformations with (independent) external measurements.
/**
 * Composed of curved trajectories in space.
 * \param [in] aPointsAndTransList List containing pairs with list of points and transformation (at inner (first) point)
 * \param [in] extDerivatives Derivatives of external measurements vs external parameters
 * \param [in] extMeasurements External measurements (residuals)
 * \param [in] extPrecisions Precision of external measurements
 */
GblTrajectory::GblTrajectory(const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> > &aPointsAndTransList,
		const TMatrixD &extDerivatives, const TVectorD &extMeasurements,
		const TVectorD &extPrecisions) :
numAllPoints(), numPoints(), numOffsetPoints(0), numOffsets(0), numInnerTransformations(aPointsAndTransList.size()),
numParameters(0), numLocals(0), numMeasurements(0), externalPoint(0),
skippedMeasLabel(0), maxNumGlobals(0), theDimension(0), thePoints(), theData(),
measDataIndex(), scatDataIndex(), externalSeed(), innerTransformations(),
externalDerivatives(), externalMeasurements(), externalPrecisions() {

	// convert from ROOT
	unsigned int nExtMeas = extDerivatives.GetNrows();
	unsigned int nExtPar = extDerivatives.GetNcols();
	externalDerivatives.resize(nExtMeas, nExtPar);
	externalMeasurements.resize(nExtMeas);
	externalPrecisions.resize(nExtMeas);
	for (unsigned int i = 0; i < nExtMeas; ++i) {
		externalMeasurements(i) = extMeasurements[i];
		externalPrecisions(i) = extPrecisions[i];
		for (unsigned int j = 0; j < nExtPar; ++j) {
			externalDerivatives(i, j) = extDerivatives(i, j);
		}
	}
	for (unsigned int iTraj = 0; iTraj < aPointsAndTransList.size(); ++iTraj) {
		thePoints.emplace_back(std::move(aPointsAndTransList[iTraj].first));
		numPoints.push_back(thePoints.back().size());
		numAllPoints += numPoints.back();
		// convert from ROOT
		unsigned int nRowTrans = aPointsAndTransList[iTraj].second.GetNrows();
		unsigned int nColTrans = aPointsAndTransList[iTraj].second.GetNcols();
		MatrixXd aTrans(nRowTrans, nColTrans);
		for (unsigned int i = 0; i < nRowTrans; ++i) {
			for (unsigned int j = 0; j < nColTrans; ++j) {
				aTrans(i, j) = aPointsAndTransList[iTraj].second(i, j);
			}
		}
		innerTransformations.emplace_back(std::move(aTrans));
	}
	theDimension.push_back(0);
	theDimension.push_back(1);
	// kinematic (2) or geometric (1) constraint
	numInnerTransOffsets = innerTransformations[0].rows() == 5 ? 2 : 1;
	numCurvature =
	innerTransformations[0].rows() == 5 ?
	innerTransformations[0].cols() :
	innerTransformations[0].cols() + numInnerTransformations;
	construct();// construct (composed) trajectory
}

/// Create new composed trajectory from list of points and transformations with (correlated) external measurements.
/**
 * Composed of curved trajectories in space.
 * \param [in] aPointsAndTransList List containing pairs with list of points and transformation (at inner (first) point)
 * \param [in] extDerivatives Derivatives of external measurements vs external parameters
 * \param [in] extMeasurements External measurements (residuals)
 * \param [in] extPrecisions Precision of external measurements
 */
GblTrajectory::GblTrajectory(const std::vector<std::pair<std::vector<GblPoint>, TMatrixD> > &aPointsAndTransList,
		const TMatrixD &extDerivatives, const TVectorD &extMeasurements,
		const TMatrixDSym &extPrecisions) :
numAllPoints(), numPoints(), numOffsetPoints(0), numOffsets(0), numInnerTransformations(aPointsAndTransList.size()),
numParameters(0), numLocals(0), numMeasurements(0), externalPoint(0),
skippedMeasLabel(0), maxNumGlobals(0), theDimension(0), thePoints(), theData(),
measDataIndex(), scatDataIndex(), externalSeed(), innerTransformations() {

	// diagonalize external measurement
	TMatrixDSymEigen extEigen(extPrecisions);
	TMatrixD extTransformation = extEigen.GetEigenVectors();
	extTransformation.T();
	TMatrixD aDerivatives = extTransformation * extDerivatives;
	TVectorD aMeasurements = extTransformation * extMeasurements;
	TVectorD aPrecisions = extEigen.GetEigenValues();
	// convert from ROOT
	unsigned int nExtMeas = aDerivatives.GetNrows();
	unsigned int nExtPar = aDerivatives.GetNcols();
	externalDerivatives.resize(nExtMeas, nExtPar);
	externalMeasurements.resize(nExtMeas);
	externalPrecisions.resize(nExtMeas);
	for (unsigned int i = 0; i < nExtMeas; ++i) {
		externalMeasurements(i) = aMeasurements[i];
		externalPrecisions(i) = aPrecisions[i];
		for (unsigned int j = 0; j < nExtPar; ++j) {
			externalDerivatives(i, j) = aDerivatives(i, j);
		}
	}
	for (unsigned int iTraj = 0; iTraj < aPointsAndTransList.size(); ++iTraj) {
		thePoints.emplace_back(std::move(aPointsAndTransList[iTraj].first));
		numPoints.push_back(thePoints.back().size());
		numAllPoints += numPoints.back();
		// convert from ROOT
		unsigned int nRowTrans = aPointsAndTransList[iTraj].second.GetNrows();
		unsigned int nColTrans = aPointsAndTransList[iTraj].second.GetNcols();
		MatrixXd aTrans(nRowTrans, nColTrans);
		for (unsigned int i = 0; i < nRowTrans; ++i) {
			for (unsigned int j = 0; j < nColTrans; ++j) {
				aTrans(i, j) = aPointsAndTransList[iTraj].second(i, j);
			}
		}
		innerTransformations.emplace_back(std::move(aTrans));
	}
	theDimension.push_back(0);
	theDimension.push_back(1);
	// kinematic (2) or geometric (1) constraint
	numInnerTransOffsets = innerTransformations[0].rows() == 5 ? 2 : 1;
	numCurvature =
	innerTransformations[0].rows() == 5 ?
	innerTransformations[0].cols() :
	innerTransformations[0].cols() + numInnerTransformations;
	construct();// construct (composed) trajectory
}
#endif

GblTrajectory::~GblTrajectory() {
}

/// Retrieve validity of trajectory
bool GblTrajectory::isValid() const {
	return constructOK;
}

/// Retrieve number of point from trajectory
unsigned int GblTrajectory::getNumPoints() const {
	return numAllPoints;
}

/// Construct trajectory from list of points.
/**
 * Trajectory is prepared for fit or output to binary file, may consists of sub-trajectories.
 */
void GblTrajectory::construct() {

	constructOK = false;
	fitOK = false;
	unsigned int aLabel = 0;
	if (numAllPoints < 2) {
		std::cout << " GblTrajectory construction failed: too few GblPoints "
				<< std::endl;
		return;
	}
	// loop over trajectories
	numTrajectories = thePoints.size();
	std::vector<GblMeasurement>::const_iterator itMeas;

	//std::cout << " numTrajectories: " << numTrajectories << ", " << innerTransformations.size() << std::endl;
	for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
		std::vector<GblPoint>::iterator itPoint;
		for (itPoint = thePoints[iTraj].begin();
				itPoint < thePoints[iTraj].end(); ++itPoint) {
			for (itMeas = itPoint->getMeasBegin();
					itMeas < itPoint->getMeasEnd(); ++itMeas) {
				if (itMeas->isEnabled()) {
					numLocals = std::max(numLocals, itMeas->getNumLocals());
					numMeasurements += itMeas->getMeasDim();
				}
			}
			itPoint->setLabel(++aLabel);
		}
	}
	defineOffsets();
	calcJacobians();
	try {
		prepare();
	} catch (std::overflow_error &e) {
		std::cout << " GblTrajectory construction failed: " << e.what()
				<< std::endl;
		return;
	} catch (int i) {
		std::cout << " GblTrajectory construction failed with exception " << i
				<< std::endl;
		return;
	}

	constructOK = true;
	// number of fit parameters
	numParameters =
			(numOffsets - numInnerTransOffsets * numInnerTransformations)
					* theDimension.size() + numCurvature + numLocals;
}

/// Define offsets from list of points.
/**
 * Define offsets at points with scatterers and first and last point.
 * All other points need interpolation from adjacent points with offsets.
 */
void GblTrajectory::defineOffsets() {

	// loop over trajectories
	for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
		// first point is offset
		thePoints[iTraj].front().setType(-1);
		thePoints[iTraj].front().setOffset(numOffsets++);
		numOffsetPoints++;
		// intermediate scatterers are offsets
		std::vector<GblPoint>::iterator itPoint;
		for (itPoint = thePoints[iTraj].begin() + 1;
				itPoint < thePoints[iTraj].end() - 1; ++itPoint) {
			unsigned int scatDim = itPoint->getScatDim();
			if (scatDim) {
				itPoint->setOffset(numOffsets++);
				numOffsetPoints++;
				// thick scatterer?
				if (scatDim == 4)
					numOffsets++;
			} else {
				itPoint->setOffset(-numOffsets);
			}
		}
		// last point is offset
		thePoints[iTraj].back().setType(1);
		thePoints[iTraj].back().setOffset(numOffsets++);
		numOffsetPoints++;
		// thick scatterer at last point?
		if (thePoints[iTraj].back().getScatDim() == 4)
			numOffsets++;
	}
}

/// Calculate Jacobians to previous/next scatterer from point to point ones.
void GblTrajectory::calcJacobians() {

	Matrix5d scatJacobian;
	// loop over trajectories
	for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
		// forward propagation (all)
		GblPoint *previousPoint = &thePoints[iTraj].front();
		unsigned int numStep = 0;
		std::vector<GblPoint>::iterator itPoint;
		for (itPoint = thePoints[iTraj].begin() + 1;
				itPoint < thePoints[iTraj].end(); ++itPoint) {
			if (numStep == 0) {
				scatJacobian = itPoint->getP2pJacobian();
			} else {
				scatJacobian = itPoint->getP2pJacobian() * scatJacobian;
			}
			numStep++;
			itPoint->addPrevJacobian(scatJacobian); // iPoint -> previous scatterer
			if (itPoint->getOffset() >= 0) {
				previousPoint->addNextJacobian(scatJacobian); // lastPoint -> next scatterer
				numStep = 0;
				previousPoint = &(*itPoint);
			}
		}
		// backward propagation (without scatterers)
		for (itPoint = thePoints[iTraj].end() - 1;
				itPoint > thePoints[iTraj].begin(); --itPoint) {
			if (itPoint->getOffset() >= 0) {
				scatJacobian = itPoint->getP2pJacobian();
				continue; // skip offsets
			}
			itPoint->addNextJacobian(scatJacobian); // iPoint -> next scatterer
			scatJacobian = scatJacobian * itPoint->getP2pJacobian();
		}
	}
}

/// Get jacobian for transformation from fit to track parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i,u_i+1..) to track (q/p,u',u) parameters
 * including additional local parameters (and transformation from external parameters).
 * \param [in] aSignedLabel (Signed) label of point for external seed
 * (<0: in front, >0: after point, slope changes at scatterer!)
 * \return List of fit parameters with non zero derivatives and
 * corresponding transformation matrix
 */
std::pair<std::vector<unsigned int>, MatrixXd> GblTrajectory::getJacobian(
		int aSignedLabel) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;
	// find trajectory, position in trajectory
	unsigned int aLabel = abs(aSignedLabel);
	unsigned int firstLabel = 1;
	unsigned int lastLabel = 0;
	unsigned int aTrajectory = 0;
	// loop over trajectories
	for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
		aTrajectory = iTraj;
		lastLabel += numPoints[iTraj];
		if (aLabel <= lastLabel)
			break;
		if (iTraj < numTrajectories - 1)
			firstLabel += numPoints[iTraj];
	}
	int nJacobian; // 0: prev, 1: next
	// check consistency of (index, direction)
	if (aSignedLabel > 0) {
		nJacobian = 1;
		if (aLabel >= lastLabel) {
			aLabel = lastLabel;
			nJacobian = 0;
		}
	} else {
		nJacobian = 0;
		if (aLabel <= firstLabel) {
			aLabel = firstLabel;
			nJacobian = 1;
		}
	}

	const GblPoint aPoint = thePoints[aTrajectory][aLabel - firstLabel];
	std::array<unsigned int, 5> labDer;
	Matrix5d matDer;
	getFitToLocalJacobian(labDer, matDer, aPoint, 5, nJacobian);

	unsigned int nBand = 0; // number of parameters from band part
	if (numInnerTransformations > 0) {
		unsigned int lastExt = innerTransLab[aTrajectory][2
				* numInnerTransOffsets]; // last label for external parameters
		for (unsigned int i = 0; i < 5; ++i)
			if (labDer[i] > lastExt)
				nBand++;
	} else {
		nBand = 2 * nDim;
	}

	unsigned int nBorder = nCurv + nLocals;
	unsigned int nParBRL = nBorder + nBand;
	unsigned int nParLoc = nLocals + 5;
	std::vector<unsigned int> anIndex;
	anIndex.reserve(nParBRL);
	MatrixXd aJacobian(nParLoc, nParBRL);
	aJacobian.setZero();

	// from local parameters
	for (unsigned int i = 0; i < nLocals; ++i) {
		aJacobian(i + 5, i) = 1.0;
		anIndex.push_back(i + 1);
	}
	// from trajectory parameters
	unsigned int iCol = nLocals;
	// composed trajectory ?
	if (numInnerTransformations > 0) {
		// transformation external to (simple) broken line (fit) parameters
		unsigned int nSimple = nCurv + nBand;
		MatrixXd aTrans(5, nSimple);
		aTrans.setZero();
		// external part, curvature
		aTrans.block(0, 0, 1, nCurv) = innerTransDer[aTrajectory].block(0, 0, 1,
				nCurv);
		for (unsigned int i = 0; i < nCurv; ++i)
			anIndex.push_back(++iCol);
		if (nBand < 4) {
			// external part, offsets (nBand=0: 2, nBand=2: 1)
			unsigned int iRow = 1 + 2 * numInnerTransOffsets + nBand - 4; // start row (1: 1st, 3: 2nd offset)
			aTrans.block(1, 0, 4 - nBand, nCurv) =
					innerTransDer[aTrajectory].block(iRow, 0, 4 - nBand, nCurv);
		}
		// (remaining) band part
		for (unsigned int i = 5 - nBand; i < 5; ++i) {
			aTrans(i, iCol++) = 1.;
			anIndex.push_back(
					labDer[i]
							- numInnerTransOffsets * nDim * (aTrajectory + 1)); // adjust label
		}
		aJacobian.block(0, nLocals, 5, nSimple) = matDer * aTrans;
	} else {
		// simple trajectory
		for (unsigned int i = 0; i < 5; ++i) {
			if (labDer[i] > 0) {
				anIndex.push_back(labDer[i]);
				for (unsigned int j = 0; j < 5; ++j) {
					aJacobian(j, iCol) = matDer(j, i);
				}
				++iCol;
			}
		}
	}
	return std::make_pair(anIndex, aJacobian);
}

/// Get (part of) jacobian for transformation from (trajectory) fit to track parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i,u_i+1..) to local (q/p,u',u) parameters.
 * \param [out] anIndex List of fit parameters (zero for zero derivatives)
 * \param [out] aJacobian Corresponding transformation matrix
 * \param [in] aPoint Point to use
 * \param [in] measDim Dimension of 'measurement'
 * (<=2: calculate only offset part, >2: complete matrix)
 * \param [in] nJacobian Direction (0: to previous offset, 1: to next offset)
 */
void GblTrajectory::getFitToLocalJacobian(std::array<unsigned int, 5> &anIndex,
		Matrix5d &aJacobian, const GblPoint &aPoint, unsigned int measDim,
		unsigned int nJacobian) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;

	int nOffset = aPoint.getOffset();
	unsigned int scatDim = aPoint.getScatDim();

	anIndex = { }; // reset to 0
	aJacobian.setZero();
	if (nOffset < 0) // need interpolation
			{
		Matrix2d prevW, prevWJ, nextW, nextWJ, matN;
		Vector2d prevWd, nextWd;
		aPoint.getDerivatives(0, prevW, prevWJ, prevWd); // W-, W- * J-, W- * d-
		aPoint.getDerivatives(1, nextW, nextWJ, nextWd); // W+, W+ * J+, W+ * d+
		const Matrix2d sumWJ(prevWJ + nextWJ);
		matN = sumWJ.inverse(); // N = (W- * J- + W+ * J+)^-1
		// derivatives for u_int
		const Matrix2d prevNW(matN * prevW); // N * W-
		const Matrix2d nextNW(matN * nextW); // N * W+
		const Vector2d prevNd(matN * prevWd); // N * W- * d-
		const Vector2d nextNd(matN * nextWd); // N * W+ * d+

		unsigned int iOff = nDim * (-nOffset - 1) + nLocals + nCurv + 1; // first offset ('i' in u_i)

		// local offset
		if (nCurv > 0) {
			aJacobian.block<2, 1>(3, 0) = -prevNd - nextNd; // from curvature
			anIndex[0] = nLocals + 1;
		}
		aJacobian.block<2, 2>(3, 1) = prevNW; // from 1st Offset
		aJacobian.block<2, 2>(3, 3) = nextNW; // from 2nd Offset
		for (unsigned int i = 0; i < nDim; ++i) {
			anIndex[1 + theDimension[i]] = iOff + i;
			anIndex[3 + theDimension[i]] = iOff + nDim + i;
		}

		// local slope and curvature
		if (measDim > 2) {
			// derivatives for u'_int
			const Matrix2d prevWPN(nextWJ * prevNW); // W+ * J+ * N * W-
			const Matrix2d nextWPN(prevWJ * nextNW); // W- * J- * N * W+
			const Vector2d prevWNd(nextWJ * prevNd); // W+ * J+ * N * W- * d-
			const Vector2d nextWNd(prevWJ * nextNd); // W- * J- * N * W+ * d+
			if (nCurv > 0) {
				aJacobian(0, 0) = 1.0;
				aJacobian.block<2, 1>(1, 0) = prevWNd - nextWNd; // from curvature
			}
			aJacobian.block<2, 2>(1, 1) = -prevWPN; // from 1st Offset
			aJacobian.block<2, 2>(1, 3) = nextWPN; // from 2nd Offset
		}
	} else { // at point
		// anIndex must be sorted
		// forward : iOff2 = iOff1 + nDim, index1 = 1, index2 = 3
		// backward: iOff2 = iOff1 - nDim, index1 = 3, index2 = 1
		// adjust for thick scatterer (before/after)
		if (scatDim == 4)
			nOffset += nJacobian;
		unsigned int iOff1 = nDim * nOffset + nCurv + nLocals + 1; // first offset ('i' in u_i)
		unsigned int index1 = 3 - 2 * nJacobian; // index of first offset
		unsigned int iOff2 = iOff1 + nDim * (nJacobian * 2 - 1); // second offset ('i' in u_i)
		unsigned int index2 = 1 + 2 * nJacobian; // index of second offset
		// local offset
		aJacobian(3, index1) = 1.0; // from 1st Offset
		aJacobian(4, index1 + 1) = 1.0;
		for (unsigned int i = 0; i < nDim; ++i) {
			anIndex[index1 + theDimension[i]] = iOff1 + i;
		}

		// local slope and curvature
		if (measDim > 2) {
			Matrix2d matW, matWJ;
			Vector2d vecWd;
			aPoint.getDerivatives(nJacobian, matW, matWJ, vecWd); // W, W * J, W * d
			double sign = (nJacobian > 0) ? 1. : -1.;
			if (nCurv > 0) {
				aJacobian(0, 0) = 1.0;
				aJacobian.block<2, 1>(1, 0) = -sign * vecWd; // from curvature
				anIndex[0] = nLocals + 1;
			}
			aJacobian.block<2, 2>(1, index1) = -sign * matWJ; // from 1st Offset
			aJacobian.block<2, 2>(1, index2) = sign * matW; // from 2nd Offset
			for (unsigned int i = 0; i < nDim; ++i) {
				anIndex[index2 + theDimension[i]] = iOff2 + i;
			}
		}
	}
}

/// Get jacobian for transformation from (trajectory) fit to kink parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i-1,u_i,u_i+1..) to kink (du') parameters.
 * \param [out] anIndex List of fit parameters (zero for zero derivatives)
 * \param [out] aJacobian Corresponding transformation matrix
 * \param [in] aPoint Point to use
 * \return Number of derivatives
 */
unsigned int GblTrajectory::getFitToKinkJacobian(
		std::array<unsigned int, 9> &anIndex, Matrix49d &aJacobian,
		const GblPoint &aPoint) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;

	int nOffset = aPoint.getOffset();

	anIndex = { }; // reset to 0
	aJacobian.topRows<2>().setZero();

	Matrix2d prevW, prevWJ, nextW, nextWJ;
	Vector2d prevWd, nextWd;
	aPoint.getDerivatives(0, prevW, prevWJ, prevWd); // W-, W- * J-, W- * d-
	aPoint.getDerivatives(1, nextW, nextWJ, nextWd); // W-, W- * J-, W- * d-
	const Matrix2d sumWJ(prevWJ + nextWJ); // W- * J- + W+ * J+
	const Vector2d sumWd(prevWd + nextWd); // W+ * d+ + W- * d-

	unsigned int iOff = (nOffset - 1) * nDim + nCurv + nLocals + 1; // first offset ('i' in u_i)

	// local offset
	if (nCurv > 0) {
		aJacobian.block<2, 1>(0, 0) = -sumWd; // from curvature
		anIndex[0] = nLocals + 1;
	}
	aJacobian.block<2, 2>(0, 1) = prevW; // from 1st Offset
	aJacobian.block<2, 2>(0, 3) = -sumWJ; // from 2nd Offset
	aJacobian.block<2, 2>(0, 5) = nextW; // from 3rd Offset
	for (unsigned int i = 0; i < nDim; ++i) {
		anIndex[1 + theDimension[i]] = iOff + i;
		anIndex[3 + theDimension[i]] = iOff + nDim + i;
		anIndex[5 + theDimension[i]] = iOff + nDim * 2 + i;
	}
	return 7;
}

/// Get jacobian for transformation from (trajectory) fit to step parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i-1,u_i-,u_i+,u_i+1..) to step (du) parameters.
 * \param [out] anIndex List of fit parameters (zero for zero derivatives)
 * \param [out] aJacobian Corresponding transformation matrix
 * \param [in] aPoint Point to use
 * \return Number of derivatives
 */
unsigned int GblTrajectory::getFitToStepJacobian(
		std::array<unsigned int, 9> &anIndex, Matrix49d &aJacobian,
		const GblPoint &aPoint) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;

	int nOffset = aPoint.getOffset();

	anIndex = { }; // reset to 0
	aJacobian.leftCols<4>().setZero();

	unsigned int iOff = (nOffset - 1) * nDim + nCurv + nLocals + 1; // first offset ('i' in u_i)

	// step
	aJacobian(2, 0) = -1.;  // from 2nd Offset
	aJacobian(3, 1) = -1.;  // from 2nd Offset
	aJacobian(2, 2) = +1.;  // from 3rd Offset
	aJacobian(3, 3) = +1.;  // from 3rd Offset

	for (unsigned int i = 0; i < nDim; ++i) {
		anIndex[0 + theDimension[i]] = iOff + nDim + i;
		anIndex[2 + theDimension[i]] = iOff + nDim * 2 + i;
	}
	return 4;
}

/// Get jacobian for transformation from (trajectory) fit to kink and step parameters at point.
/**
 * Jacobian broken lines (q/p,..,u_i-1,u_i-,u_i+,u_i+1..) to kink (du') and step (du) parameters.
 * \param [out] anIndex List of fit parameters (zero for zero derivatives)
 * \param [out] aJacobian Corresponding transformation matrix
 * \param [in] aPoint Point to use
 * \return Number of derivatives
 */
unsigned int GblTrajectory::getFitToKinkAndStepJacobian(
		std::array<unsigned int, 9> &anIndex, Matrix49d &aJacobian,
		const GblPoint &aPoint) const {

	unsigned int nDim = theDimension.size();
	unsigned int nCurv = numCurvature;
	unsigned int nLocals = numLocals;

	int nOffset = aPoint.getOffset();

	anIndex = { }; // reset to 0
	aJacobian.setZero();

	Matrix2d prevW, prevWJ, nextW, nextWJ;
	Vector2d prevWd, nextWd;
	aPoint.getDerivatives(0, prevW, prevWJ, prevWd); // W-, W- * J-, W- * d-
	aPoint.getDerivatives(1, nextW, nextWJ, nextWd); // W-, W- * J-, W- * d-

	unsigned int iOff = (nOffset - 1) * nDim + nCurv + nLocals + 1; // first offset ('i' in u_i)

	// kink
	if (nCurv > 0) {
		aJacobian.block<2, 1>(0, 0) = -(prevWd + nextWd); // from curvature
		anIndex[0] = nLocals + 1;
	}
	aJacobian.block<2, 2>(0, 1) = prevW; // from 1st Offset
	aJacobian.block<2, 2>(0, 3) = -prevWJ; // from 2nd Offset
	aJacobian.block<2, 2>(0, 5) = -nextWJ; // from 3rd Offset
	aJacobian.block<2, 2>(0, 7) = nextW; // from 4th Offset
	// step
	aJacobian(2, 3) = -1.;  // from 2nd Offset
	aJacobian(3, 4) = -1.;  // from 2nd Offset
	aJacobian(2, 5) = +1.;  // from 3rd Offset
	aJacobian(3, 6) = +1.;  // from 3rd Offset

	for (unsigned int i = 0; i < nDim; ++i) {
		anIndex[1 + theDimension[i]] = iOff + i;
		anIndex[3 + theDimension[i]] = iOff + nDim + i;
		anIndex[5 + theDimension[i]] = iOff + nDim * 2 + i;
		anIndex[7 + theDimension[i]] = iOff + nDim * 3 + i;
	}
	return 9;
}

/// Get fit results for external parameters.
/**
 * Get corrections and covariance matrix for external parameters.
 *
 * \param [out] extPar Corrections for external parameters
 * \param [out] extCov Covariance for external parameters
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getExtResults(Eigen::VectorXd &extPar,
		Eigen::MatrixXd &extCov) const {
	if (not fitOK)
		return 1;
	// get external parameters (single block after locals)
	unsigned int nExt = innerTransformations[0].cols();
	VectorXd aVec(nExt); // compressed vector
	std::vector<unsigned int> index;
	for (unsigned int i = 0; i < nExt; ++i) {
		aVec[i] = theVector(numLocals + i);
		index.push_back(numLocals + i + 1);
	}
	extPar = aVec;
	extCov = theMatrix.getBlockMatrix(index); // compressed matrix
	return 0;
}

/// Get fit results at point.
/**
 * Get corrections and covariance matrix for local track and additional parameters
 * in forward or backward direction.
 *
 * The point is identified by its label (1..number(points)), the sign distinguishes the
 * backward (facing previous point) and forward 'side' (facing next point).
 * For (thick) scatterers the track direction (and offset) may change in between.
 *
 * \param [in] aSignedLabel (Signed) label of point on trajectory
 * (<0: in front, >0: after point, slope (and offset) changes at (thick) scatterer!)
 * \param [out] localPar Corrections for local parameters
 * \param [out] localCov Covariance for local parameters
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getResults(int aSignedLabel,
		Eigen::VectorXd &localPar, Eigen::MatrixXd &localCov) const {
	if (not fitOK)
		return 1;
	std::pair<std::vector<unsigned int>, MatrixXd> indexAndJacobian =
			getJacobian(aSignedLabel);
	unsigned int nParBrl = indexAndJacobian.first.size();
	VectorXd aVec(nParBrl); // compressed vector
	for (unsigned int i = 0; i < nParBrl; ++i) {
		aVec[i] = theVector(indexAndJacobian.first[i] - 1);
	}
	MatrixXd aMat = theMatrix.getBlockMatrix(indexAndJacobian.first); // compressed matrix
	localPar = indexAndJacobian.second * aVec;
	localCov = indexAndJacobian.second * aMat
			* indexAndJacobian.second.transpose();
	return 0;
}

/// Get residuals from fit at point for measurement "long list".
/**
 * Get (diagonalized) residual, error of measurement and residual and down-weighting
 * factor for measurement at point
 *
 * \param [in]  aLabel Label of point on trajectory
 * \param [out] numData Number of data blocks from measurement at point
 * \param [out] aResiduals Measurements-Predictions
 * \param [out] aMeasErrors Errors of Measurements
 * \param [out] aResErrors Errors of Residuals (including correlations from track fit)
 * \param [out] aDownWeights Down-Weighting factors
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getMeasResults(unsigned int aLabel,
		unsigned int &numData, Eigen::VectorXd &aResiduals,
		Eigen::VectorXd &aMeasErrors, Eigen::VectorXd &aResErrors,
		Eigen::VectorXd &aDownWeights) {
	numData = 0;
	if (not fitOK)
		return 1;

	unsigned int firstData = measDataIndex[aLabel - 1]; // first data block with measurement
	numData = measDataIndex[aLabel] - firstData; // number of data blocks
	for (unsigned int i = 0; i < numData; ++i) {
		getResAndErr(firstData + i, (aLabel != skippedMeasLabel), aResiduals(i),
				aMeasErrors(i), aResErrors(i), aDownWeights(i));
	}
	return 0;
}

/// Get residuals from fit at point for measurement "short list".
/**
 * Get (diagonalized) residual, error of measurement and residual and down-weighting
 * factor for measurement at point
 *
 * \param [in]  aLabel Label of point on trajectory
 * \param [out] numData Number of data blocks from measurement at point
 * \param [out] aResiduals Measurements-Predictions
 * \param [out] aMeasErrors Errors of Measurements
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getMeasResults(unsigned int aLabel,
		unsigned int &numData, Eigen::VectorXd &aResiduals,
		Eigen::VectorXd &aMeasErrors) {
	numData = 0;
	if (not fitOK)
		return 1;

	unsigned int firstData = measDataIndex[aLabel - 1]; // first data block with measurement
	numData = measDataIndex[aLabel] - firstData; // number of data blocks
	for (unsigned int i = 0; i < numData; ++i) {
		getResAndErr(firstData + i, aResiduals(i), aMeasErrors(i));
	}
	return 0;
}

/// Get (kink) residuals from fit at point for scatterer.
/**
 * Get (diagonalized) residual, error of measurement and residual and down-weighting
 * factor for scatterering kinks at point
 *
 * \param [in]  aLabel Label of point on trajectory
 * \param [out] numData Number of data blocks from scatterer at point
 * \param [out] aResiduals (kink)Measurements-(kink)Predictions
 * \param [out] aMeasErrors Errors of (kink)Measurements
 * \param [out] aResErrors Errors of Residuals (including correlations from track fit)
 * \param [out] aDownWeights Down-Weighting factors
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getScatResults(unsigned int aLabel,
		unsigned int &numData, Eigen::VectorXd &aResiduals,
		Eigen::VectorXd &aMeasErrors, Eigen::VectorXd &aResErrors,
		Eigen::VectorXd &aDownWeights) {
	numData = 0;
	if (not fitOK)
		return 1;

	unsigned int firstData = scatDataIndex[aLabel - 1]; // first data block with scatterer
	numData = scatDataIndex[aLabel] - firstData; // number of data blocks
	for (unsigned int i = 0; i < numData; ++i) {
		getResAndErr(firstData + i, true, aResiduals(i), aMeasErrors(i),
				aResErrors(i), aDownWeights(i));
	}
	return 0;
}

#ifdef GBL_EIGEN_SUPPORT_ROOT
/// Get fit results for external parameters.
/**
 * Get corrections and covariance matrix for external parameters.
 *
 * \param [out] extPar Corrections for external parameters
 * \param [out] extCov Covariance for external parameters
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getExtResults(TVectorD &extPar,
		TMatrixDSym &extCov) const {
	if (not fitOK)
		return 1;
	// get external parameters (single block after locals)
	unsigned int nExt = innerTransformations[0].cols();
	VectorXd aVec(nExt); // compressed vector
	std::vector<unsigned int> index;
	for (unsigned int i = 0; i < nExt; ++i) {
		aVec[i] = theVector(numLocals + i);
		index.push_back(numLocals + i + 1);
	}
	MatrixXd aMat = theMatrix.getBlockMatrix(index); // compressed matrix
	// convert to ROOT
	unsigned int nParOut = extPar.GetNrows();
	for (unsigned int i = 0; i < nParOut; ++i) {
		extPar[i] = aVec(i);
		for (unsigned int j = 0; j < nParOut; ++j) {
			extCov(i, j) = aMat(i, j);
		}
	}
	return 0;
}

/// Get fit results at point.
/**
 * Get corrections and covariance matrix for local track and additional parameters
 * in forward or backward direction.
 *
 * The point is identified by its label (1..number(points)), the sign distinguishes the
 * backward (facing previous point) and forward 'side' (facing next point).
 * For (thick) scatterers the track direction (and offset) may change in between.
 *
 * \param [in] aSignedLabel (Signed) label of point on trajectory
 * (<0: in front, >0: after point, slope changes at scatterer!)
 * \param [out] localPar Corrections for local parameters
 * \param [out] localCov Covariance for local parameters
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getResults(int aSignedLabel, TVectorD &localPar,
		TMatrixDSym &localCov) const {
	if (not fitOK)
	return 1;
	std::pair<std::vector<unsigned int>, MatrixXd> indexAndJacobian =
	getJacobian(aSignedLabel);
	unsigned int nParBrl = indexAndJacobian.first.size();
	VectorXd aVec(nParBrl); // compressed vector
	for (unsigned int i = 0; i < nParBrl; ++i) {
		aVec[i] = theVector(indexAndJacobian.first[i] - 1);
	}
	MatrixXd aMat = theMatrix.getBlockMatrix(indexAndJacobian.first); // compressed matrix
	VectorXd aLocalPar = indexAndJacobian.second * aVec;
	MatrixXd aLocalCov = indexAndJacobian.second * aMat
	* indexAndJacobian.second.transpose();
	// convert to ROOT
	unsigned int nParOut = localPar.GetNrows();
	for (unsigned int i = 0; i < nParOut; ++i) {
		localPar[i] = aLocalPar(i);
		for (unsigned int j = 0; j < nParOut; ++j) {
			localCov(i, j) = aLocalCov(i, j);
		}
	}
	return 0;
}

/// Get residuals from fit at point for measurement.
/**
 * Get (diagonalized) residual, error of measurement and residual and down-weighting
 * factor for measurement at point
 *
 * \param [in]  aLabel Label of point on trajectory
 * \param [out] numData Number of data blocks from measurement at point
 * \param [out] aResiduals Measurements-Predictions
 * \param [out] aMeasErrors Errors of Measurements
 * \param [out] aResErrors Errors of Residuals (including correlations from track fit)
 * \param [out] aDownWeights Down-Weighting factors
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getMeasResults(unsigned int aLabel,
		unsigned int &numData, TVectorD &aResiduals, TVectorD &aMeasErrors,
		TVectorD &aResErrors, TVectorD &aDownWeights) {
	numData = 0;
	if (not fitOK)
	return 1;

	unsigned int firstData = measDataIndex[aLabel - 1]; // first data block with measurement
	numData = measDataIndex[aLabel] - firstData;// number of data blocks
	for (unsigned int i = 0; i < numData; ++i) {
		getResAndErr(firstData + i, (aLabel != skippedMeasLabel), aResiduals[i],
				aMeasErrors[i], aResErrors[i], aDownWeights[i]);
	}
	return 0;
}

/// Get (kink) residuals from fit at point for scatterer.
/**
 * Get (diagonalized) residual, error of measurement and residual and down-weighting
 * factor for scatterering kinks at point
 *
 * \param [in]  aLabel Label of point on trajectory
 * \param [out] numData Number of data blocks from scatterer at point
 * \param [out] aResiduals (kink)Measurements-(kink)Predictions
 * \param [out] aMeasErrors Errors of (kink)Measurements
 * \param [out] aResErrors Errors of Residuals (including correlations from track fit)
 * \param [out] aDownWeights Down-Weighting factors
 * \return error code (non-zero if trajectory not fitted successfully)
 */
unsigned int GblTrajectory::getScatResults(unsigned int aLabel,
		unsigned int &numData, TVectorD &aResiduals, TVectorD &aMeasErrors,
		TVectorD &aResErrors, TVectorD &aDownWeights) {
	numData = 0;
	if (not fitOK)
	return 1;

	unsigned int firstData = scatDataIndex[aLabel - 1]; // first data block with scatterer
	numData = scatDataIndex[aLabel] - firstData;// number of data blocks
	for (unsigned int i = 0; i < numData; ++i) {
		getResAndErr(firstData + i, true, aResiduals[i], aMeasErrors[i],
				aResErrors[i], aDownWeights[i]);
	}
	return 0;
}

#endif

/// Get (list of) labels of points on (simple) valid trajectory
/**
 * \param [out] aLabelList List of labels (aLabelList[i] = i+1)
 * \return error code (non-zero if trajectory not valid (constructed successfully))
 */
unsigned int GblTrajectory::getLabels(
		std::vector<unsigned int> &aLabelList) const {
	if (not constructOK)
		return 1;

	unsigned int aLabel = 0;
	unsigned int nPoint = thePoints[0].size();
	aLabelList.resize(nPoint);
	for (unsigned i = 0; i < nPoint; ++i) {
		aLabelList[i] = ++aLabel;
	}
	return 0;
}

/// Get (list of lists of) labels of points on (composed) valid trajectory
/**
 * \param [out] aLabelList List of of lists of labels
 * \return error code (non-zero if trajectory not valid (constructed successfully))
 */
unsigned int GblTrajectory::getLabels(
		std::vector<std::vector<unsigned int> > &aLabelList) const {
	if (not constructOK)
		return 1;

	unsigned int aLabel = 0;
	aLabelList.resize(numTrajectories);
	for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
		unsigned int nPoint = thePoints[iTraj].size();
		aLabelList[iTraj].resize(nPoint);
		for (unsigned i = 0; i < nPoint; ++i) {
			aLabelList[iTraj][i] = ++aLabel;
		}
	}
	return 0;
}

/// Get residual and errors from data block "long list".
/**
 * Get residual, error of measurement and residual and down-weighting
 * factor for (single) data block
 * \param [in]  aData Label of data block
 * \param [in]  used  Flag for usage of data block in fit
 * \param [out] aResidual Measurement-Prediction
 * \param [out] aMeasError Error of Measurement
 * \param [out] aResError Error of Residual (including correlations from track fit)
 * \param [out] aDownWeight Down-Weighting factor
 */
void GblTrajectory::getResAndErr(unsigned int aData, bool used,
		double &aResidual, double &aMeasError, double &aResError,
		double &aDownWeight) {

	double aMeasVar;
	unsigned int numLocal;
	unsigned int *indLocal;
	double *derLocal;
	theData[aData].getResidual(aResidual, aMeasVar, aDownWeight, numLocal,
			indLocal, derLocal);
	VectorXd aVec(numLocal); // compressed vector of derivatives
	for (unsigned int j = 0; j < numLocal; ++j) {
		aVec[j] = derLocal[j];
	}
	MatrixXd aMat = theMatrix.getBlockMatrix(numLocal, indLocal); // compressed (covariance) matrix
	double aFitVar = aVec.transpose() * aMat * aVec; // variance from track fit
	aFitVar *= aDownWeight; // account for down-weighting (of measurement in fit)
	aMeasError = sqrt(aMeasVar); // error of measurement
	if (used)
		aResError = (aFitVar < aMeasVar ? sqrt(aMeasVar - aFitVar) : 0.); // error of (biased) residual
	else
		aResError = sqrt(aMeasVar + aFitVar); // error of (unbiased) residual
}

/// Get residual and errors from data block "short list".
/**
 * Get residual, error of measurement and residual and down-weighting
 * factor for (single) data block
 * \param [in]  aData Label of data block
 * \param [out] aResidual Measurement-Prediction
 * \param [out] aMeasError Error of Measurement
 */
void GblTrajectory::getResAndErr(unsigned int aData, double &aResidual,
		double &aMeasError) {

	double aMeasVar;
	theData[aData].getResidual(aResidual, aMeasVar);
	aMeasError = sqrt(aMeasVar); // error of measurement
}

/// Build linear equation system from data (blocks).
void GblTrajectory::buildLinearEquationSystem() {
	unsigned int nBorder = numCurvature + numLocals;
	theVector.resize(numParameters);
	theVector.setZero();
	theMatrix.resize(numParameters, nBorder);
	theMatrix.setZero();
	double aValue, aWeight;
	unsigned int *indLocal;
	double *derLocal;
	unsigned int numLocal;

	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); ++itData) {
		// skipped (internal) measurement ?
		if (itData->getLabel() == skippedMeasLabel
				&& itData->getType() == InternalMeasurement)
			continue;

		itData->getLocalData(aValue, aWeight, numLocal, indLocal, derLocal);
		for (unsigned int j = 0; j < numLocal; ++j) {
			theVector(indLocal[j] - 1) += derLocal[j] * aWeight * aValue;
		}
		theMatrix.addBlockMatrix(aWeight, numLocal, indLocal, derLocal);
	}
}

/// Prepare fit for simple or composed trajectory
/**
 * Generate data (blocks) from measurements, kinks, external seed and measurements.
 *
 * \exception 10 : inner transformation matrix with invalid number of rows (valid are 5=kinematic or 2=geometric constraint)
 * \exception 11 : inner transformation matrix with too few columns (must be >= number of rows)
 * \exception 12 : inner transformation matrices with varying sizes
 * \exception 13 : too many external derivatives (must be <= number of columns of inner transformation matrix)
 */
void GblTrajectory::prepare() {
	unsigned int nDim = theDimension.size();
	// upper limit
	unsigned int maxData = numMeasurements + nDim * (numOffsets - 2)
			+ externalSeed.rows();
	theData.reserve(maxData);
	measDataIndex.resize(numAllPoints + 3); // include external seed and measurements
	scatDataIndex.resize(numAllPoints + 1);
	unsigned int nData = 0;
	// composed trajectory ?
	if (numInnerTransformations > 0) {
		unsigned int nInnerTransRows = innerTransformations[0].rows();
		unsigned int nInnerTransCols = innerTransformations[0].cols();
		// std::cout << "composed trajectory, inner transformation (" << nInnerTransRows << "," << nInnerTransCols << ")" << std::endl;
		// check size
		if (nInnerTransRows != 5 and nInnerTransRows != 2) {
			std::cout
					<< " GblTrajectory::prepare composed trajectory with bad inner transformation matrix("
					<< nInnerTransRows << "," << nInnerTransCols
					<< "): invalid number of rows" << std::endl;
			throw 10;
		}
		if (nInnerTransRows > nInnerTransCols) {
			std::cout
					<< " GblTrajectory::prepare composed trajectory with bad inner transformation matrix("
					<< nInnerTransRows << "," << nInnerTransCols
					<< "): too few columns" << std::endl;
			throw 11;
		}
		for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
			// check size
			if (nInnerTransRows != innerTransformations[iTraj].rows()
					or nInnerTransCols != innerTransformations[iTraj].cols()) {
				std::cout
						<< " GblTrajectory::prepare composed trajectory with bad inner transformation matrix["
						<< iTraj << "]: different size as [0]" << std::endl;
				throw 12;
			}
			// innermost point
			GblPoint *innerPoint = &thePoints[iTraj].front();
			// transformation fit to local track parameters
			std::array<unsigned int, 5> firstLabels;
			Matrix5d matFitToLocal;
			getFitToLocalJacobian(firstLabels, matFitToLocal, *innerPoint, 5);
			if (nInnerTransRows == 5) {
				// (full) kinematic constraint
				// transformation local track to fit parameters
				Matrix5d matLocalToFit = matFitToLocal.inverse();
				// transformation external to fit parameters at inner (first) point
				innerTransDer.emplace_back(
						matLocalToFit * innerTransformations[iTraj]);
			} else {
				// geometric constraint (including individual curvature correction per trajectory)
				MatrixXd matInnerToFit(5, numCurvature);
				matInnerToFit.setZero();
				matInnerToFit(0, nInnerTransCols + iTraj) = 1.; // curvature correction for 'iTraj'
				matInnerToFit.block(1, 0, 2, nInnerTransCols) =
						innerTransformations[iTraj]; // geometric derivatives
				innerTransDer.emplace_back(matInnerToFit);
			}
			innerTransLab.push_back(firstLabels);
		}
	}

	Matrix5d matP;              // measurements
	std::vector<GblPoint>::iterator itPoint;
	std::vector<GblMeasurement>::const_iterator itMeas;
	// limit the scope of proDer:
	{
		// transform for external parameters
		Eigen::Matrix<double, Eigen::Dynamic, 5, Eigen::ColMajor /* default */,
				5, 5> proDer;
		// loop over trajectories
		for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
			for (itPoint = thePoints[iTraj].begin();
					itPoint < thePoints[iTraj].end(); ++itPoint) {
				Vector5d aMeas, aPrec;
				unsigned int nLabel = itPoint->getLabel();
				if (itPoint->numMeasurements()) {
					for (itMeas = itPoint->getMeasBegin();
							itMeas < itPoint->getMeasEnd(); ++itMeas) {
						// skip disabled measurements
						if (!itMeas->isEnabled())
							continue;
						unsigned int measDim = itMeas->getMeasDim();
						const MatrixXd localDer = itMeas->getLocalDerivatives();
						maxNumGlobals = std::max(maxNumGlobals,
								itMeas->getNumGlobals());
						MatrixXd transDer;
						itMeas->getMeasurement(matP, aMeas, aPrec);
						double minPrecision = itMeas->getMeasPrecMin();
						unsigned int iOff = 5 - measDim; // first active component
						std::array<unsigned int, 5> labDer;
						Matrix5d matDer, matPDer;
						matPDer.topRows(iOff).setZero(); // clear unused part
						if (measDim > 2) {
							unsigned int nJacobian =
									(itPoint->isLast()) ? 0 : 1; // last point needs backward propagation (for slopes)
							getFitToLocalJacobian(labDer, matDer, *itPoint,
									measDim, nJacobian);
							matPDer = matP * matDer;
							matPDer.bottomRows(measDim) =
									matP.bottomRightCorner(measDim, measDim)
											* matDer.bottomRows(measDim);
						} else { // 'shortcut' for position measurements
							getFitToLocalJacobian(labDer, matDer, *itPoint,
									measDim, 1); // forward propagation (-> after THICK scatterer)
							matPDer.bottomRows<2>() = matP.bottomRightCorner<2,
									2>() * matDer.bottomRows<2>();
						}

						if (numInnerTransformations > 0) {
							// transform for external parameters
							proDer.resize(measDim, Eigen::NoChange);
							proDer.setZero();
							// match parameters
							unsigned int ifirst = 0;
							unsigned int ilast = 2 * numInnerTransOffsets;
							unsigned int ilabel = 0;
							unsigned int numRelated = 0;
							while (ilabel < 5) {
								if (labDer[ilabel] > 0) {
									while (ifirst <= ilast
											and innerTransLab[iTraj][ifirst]
													!= labDer[ilabel]) {
										++ifirst;
									}
									if (ifirst > ilast) {
										labDer[ilabel] -= numInnerTransOffsets
												* nDim * (iTraj + 1); // adjust label
									} else {
										// match
										labDer[ilabel] = 0; // mark as related to external parameters
										numRelated++;
										for (unsigned int k = iOff; k < 5;
												++k) {
											proDer(k - iOff, ifirst) = matPDer(
													k, ilabel);
										}
									}
								}
								++ilabel;
							}
							if (numRelated > 0) {
								transDer.resize(measDim, numCurvature);
								transDer = proDer * innerTransDer[iTraj];
							}
						}
						for (unsigned int i = iOff; i < 5; ++i) {
							if (aPrec(i) > minPrecision) {
								GblData aData(nLabel, InternalMeasurement,
										aMeas(i), aPrec(i), iTraj,
										itPoint - thePoints[iTraj].begin(),
										itMeas - itPoint->getMeasBegin());
								aData.addDerivatives(i, labDer, matPDer, iOff,
										localDer, numLocals, transDer);
								theData.emplace_back(std::move(aData));
								nData++;
							}
						}

					}
				}
				measDataIndex[nLabel] = nData;
			}
		}
	} // end of scope for proDer

	Matrix4d matT;              // kinks
	// limit the scope of proDer:
	{
		// transform for external parameters
		Eigen::Matrix<double, Eigen::Dynamic, 5, Eigen::ColMajor /* default */,
				5, 5> proDer;
		scatDataIndex[0] = nData;
		scatDataIndex[1] = nData;
		// loop over trajectories
		for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
			for (itPoint = thePoints[iTraj].begin() + 1;
					itPoint < thePoints[iTraj].end(); ++itPoint) {
				Vector4d aMeas, aPrec;
				unsigned int nLabel = itPoint->getLabel();
				unsigned int scatDim = itPoint->getScatDim();
				if (scatDim) {
					MatrixXd transDer;
					std::array<unsigned int, 9> labDer;
					Matrix49d matDer, matTDer;
					unsigned int numDer;
					if (scatDim == 4) {
						// thick scatterer, last point?
						if (itPoint->isLast()) {
							// steps (only)
							itPoint->getReducedScatterer(matT, aMeas, aPrec);
							numDer = getFitToStepJacobian(labDer, matDer,
									*itPoint);
							matTDer.leftCols<4>() = matT * matDer.leftCols<4>(); // numDer == 4 !
						} else {
							// kinks+steps
							itPoint->getScatterer(matT, aMeas, aPrec);
							numDer = getFitToKinkAndStepJacobian(labDer, matDer,
									*itPoint);
							matTDer = matT * matDer; // numDer == 9 !
						}
					} else {
						// thin scatterer, last point?
						if (itPoint->isLast())
							break;
						// kinks
						itPoint->getScatterer(matT, aMeas, aPrec);
						numDer = getFitToKinkJacobian(labDer, matDer, *itPoint);
						matTDer.topLeftCorner<2, 7>() =
								matT.topLeftCorner<2, 2>()
										* matDer.topLeftCorner<2, 7>(); // numDer == 7 !
					}
					if (numInnerTransformations > 0) {
						// transform for external parameters
						proDer.resize(scatDim, Eigen::NoChange);
						proDer.setZero();
						// match parameters
						unsigned int ifirst = 0;
						unsigned int ilast = 2 * numInnerTransOffsets;
						unsigned int ilabel = 0;
						unsigned int numRelated = 0;
						while (ilabel < numDer) {
							if (labDer[ilabel] > 0) {
								while (innerTransLab[iTraj][ifirst]
										!= labDer[ilabel] and ifirst <= ilast) {
									++ifirst;
								}
								if (ifirst > ilast) {
									labDer[ilabel] -= numInnerTransOffsets
											* nDim * (iTraj + 1); // adjust label
								} else {
									// match
									labDer[ilabel] = 0; // mark as related to external parameters
									numRelated++;
									for (unsigned int k = 0; k < scatDim; ++k) {
										proDer(k, ifirst) = matTDer(k, ilabel);
									}
								}
							}
							++ilabel;
						}
						if (numRelated > 0) {
							transDer.resize(scatDim, numCurvature);
							transDer = proDer * innerTransDer[iTraj];
						}
					}
					// loop over kinks and steps (if any)
					for (unsigned int i = 0; i < scatDim; ++i) {
						if (aPrec(i) > 0.) {
							GblData aData(nLabel, InternalKink, aMeas(i),
									aPrec(i), iTraj,
									itPoint - thePoints[iTraj].begin());
							aData.addDerivatives(i, numDer, labDer, matTDer,
									numLocals, transDer);
							theData.emplace_back(std::move(aData));
							nData++;
						}
					}
				}
				scatDataIndex[nLabel] = nData;
			}
			scatDataIndex[thePoints[iTraj].back().getLabel()] = nData;
		}
	}

	// external seed
	if (externalPoint) {
		std::pair<std::vector<unsigned int>, MatrixXd> indexAndJacobian =
				getJacobian(externalPoint);
		std::vector<unsigned int> externalSeedIndex = indexAndJacobian.first;
		std::vector<double> externalSeedDerivatives(externalSeedIndex.size());
		SelfAdjointEigenSolver<MatrixXd> externalSeedEigen(externalSeed);
		VectorXd valEigen = externalSeedEigen.eigenvalues();
		MatrixXd vecEigen = externalSeedEigen.eigenvectors();
		vecEigen = vecEigen.transpose() * indexAndJacobian.second;
		for (int i = 0; i < externalSeed.rows(); ++i) {
			if (valEigen(i) > 0.) {
				for (int j = 0; j < externalSeed.cols(); ++j) {
					externalSeedDerivatives[j] = vecEigen(i, j);
				}
				GblData aData(externalPoint, ExternalSeed, 0., valEigen(i));
				aData.addDerivatives(externalSeedIndex,
						externalSeedDerivatives);
				theData.emplace_back(std::move(aData));
				nData++;
			}
		}
	}
	measDataIndex[numAllPoints + 1] = nData;
	// external measurements
	unsigned int nExt = externalMeasurements.rows();
	if (nExt > 0) {
		unsigned int nInnerTransCols = innerTransformations[0].cols();
		unsigned int nExtDer = externalDerivatives.cols();
		if (nExtDer > nInnerTransCols) {
			std::cout
					<< " GblTrajectory::prepare external measurement with too many derivatives: "
					<< nExtDer << ", defined: " << nInnerTransCols << std::endl;
			throw 13;
		}
		std::vector<unsigned int> index(nExtDer);
		std::vector<double> derivatives(nExtDer);
		for (unsigned int iExt = 0; iExt < nExt; ++iExt) {
			for (unsigned int iCol = 0; iCol < nExtDer; ++iCol) {
				index[iCol] = numLocals + iCol + 1;
				derivatives[iCol] = externalDerivatives(iExt, iCol);
			}
			GblData aData(1U, ExternalMeasurement, externalMeasurements(iExt),
					externalPrecisions(iExt));
			aData.addDerivatives(index, derivatives);
			theData.emplace_back(std::move(aData));
			nData++;
		}
	}
	measDataIndex[numAllPoints + 2] = nData;
}

/// Calculate predictions for all points.
void GblTrajectory::predict() {
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); ++itData) {
		itData->setPrediction(theVector);
	}
}

/// Down-weight all points.
/**
 * \param [in] aMethod M-estimator (1: Tukey, 2:Huber, 3:Cauchy)
 */
double GblTrajectory::downWeight(unsigned int aMethod) {
	double aLoss = 0.;
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData < theData.end(); ++itData) {
		aLoss += (1. - itData->setDownWeighting(aMethod));
	}
	return aLoss;
}

/// Perform fit of (valid) trajectory.
/**
 * Optionally iterate for outlier down-weighting.
 * Fit may fail due to singular or not positive definite matrices (internal exceptions 1-3).
 *
 * \param [out] Chi2 Chi2 sum (corrected for down-weighting)
 * \param [out] Ndf  Number of degrees of freedom
 * \param [out] lostWeight Sum of weights lost due to down-weighting
 * \param [in] optionList Iterations for down-weighting
 * (One character per iteration: t,h,c (or T,H,C) for Tukey, Huber or Cauchy function)
 * \param [in] aLabel Label of point where to skip measurements (for unbiased residuals)
 * \return Error code (non zero value indicates failure of fit)
 */
unsigned int GblTrajectory::fit(double &Chi2, int &Ndf, double &lostWeight,
		const std::string &optionList, unsigned int aLabel) {
	const double normChi2[4] = { 1.0, 0.8737, 0.9326, 0.8228 };
	const std::string methodList = "TtHhCc";

	Chi2 = 0.;
	Ndf = -1;
	lostWeight = 0.;
	if (not constructOK)
		return 10;

	unsigned int aMethod = 0;
	skippedMeasLabel = aLabel;

	buildLinearEquationSystem();
	lostWeight = 0.;
	unsigned int ierr = 0;
	try {

		theMatrix.solveAndInvertBorderedBand(theVector, theVector);
		predict();

		for (unsigned int i = 0; i < optionList.size(); ++i) // down weighting iterations
				{
			size_t aPosition = methodList.find(optionList[i]);
			if (aPosition != std::string::npos) {
				aMethod = aPosition / 2 + 1;
				lostWeight = downWeight(aMethod);
				buildLinearEquationSystem();
				theMatrix.solveAndInvertBorderedBand(theVector, theVector);
				predict();
			}
		}
		Ndf = -numParameters;
		Chi2 = 0.;
		for (unsigned int i = 0; i < theData.size(); ++i) {
			// skipped (internal) measurement ?
			if (theData[i].getLabel() == skippedMeasLabel
					&& theData[i].getType() == InternalMeasurement)
				continue;
			Chi2 += theData[i].getChi2();
			Ndf++;
		}
		Chi2 /= normChi2[aMethod];
		fitOK = true;

	} catch (int e) {
		std::cout << " GblTrajectory::fit exception " << e << std::endl;
		ierr = e;
	}
	return ierr;
}

/// Write valid trajectory to Millepede-II binary file.
/**
 * Trajectory state after construction (independent of fitting) is used.
 */
void GblTrajectory::milleOut(MilleBinary &aMille) {
	double aValue;
	double aErr;
	unsigned int aTraj;
	unsigned int aPoint;
	unsigned int aMeas;
	unsigned int aRow;
	unsigned int numLocal;
	unsigned int *labLocal;
	double *derLocal;
	std::vector<int> labGlobal;
	std::vector<double> derGlobal;

	if (not constructOK)
		return;

	//   data: measurements, kinks and external seed
	labGlobal.reserve(maxNumGlobals);
	derGlobal.reserve(maxNumGlobals);
	std::vector<GblData>::iterator itData;
	for (itData = theData.begin(); itData != theData.end(); ++itData) {
		itData->getAllData(aValue, aErr, numLocal, labLocal, derLocal, aTraj,
				aPoint, aMeas, aRow);
		if (itData->getType() == InternalMeasurement)
			thePoints[aTraj][aPoint].getGlobalLabelsAndDerivatives(aMeas, aRow,
					labGlobal, derGlobal);
		else
			labGlobal.resize(0);
		aMille.addData(aValue, aErr, numLocal, labLocal, derLocal, labGlobal,
				derGlobal);
	}
	aMille.writeRecord();
}

/// Print GblTrajectory
/**
 * \param [in] level print level (0: minimum, >0: more)
 */
void GblTrajectory::printTrajectory(unsigned int level) const {
	if (numInnerTransformations) {
		std::cout << "Composed GblTrajectory, " << numInnerTransformations
				<< " subtrajectories, type " << numInnerTransOffsets
				<< std::endl;
	} else {
		std::cout << "Simple GblTrajectory" << std::endl;
	}
	if (theDimension.size() < 2) {
		std::cout << " 2D-trajectory" << std::endl;
	}
	std::cout << " Number of GblPoints          : " << numAllPoints
			<< std::endl;
	std::cout << " Number of points with offsets: " << numOffsetPoints
			<< std::endl;
	std::cout << " Number of (1D or 2D) offsets : " << numOffsets << std::endl;
	std::cout << " Number of fit parameters     : " << numParameters
			<< std::endl;
	std::cout << " Number of measurements       : " << numMeasurements
			<< std::endl;
	if (externalMeasurements.rows()) {
		std::cout << " Number of ext. measurements  : "
				<< externalMeasurements.rows() << std::endl;
	}
	if (externalPoint) {
		std::cout << " Label of point with ext. seed: " << externalPoint
				<< std::endl;
	}
	if (constructOK) {
		std::cout << " Constructed OK " << std::endl;
	}
	if (fitOK) {
		std::cout << " Fitted OK " << std::endl;
	}
	if (level > 0) {
		IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
		if (numInnerTransformations) {
			std::cout << " Inner transformations" << std::endl;
			for (unsigned int i = 0; i < numInnerTransformations; ++i) {
				std::cout << innerTransformations[i].format(CleanFmt)
						<< std::endl;
			}
		}
		if (externalMeasurements.rows()) {
			std::cout << " External measurements" << std::endl;
			std::cout << "  Measurements:" << std::endl;
			std::cout << externalMeasurements.format(CleanFmt) << std::endl;
			std::cout << "  Precisions:" << std::endl;
			std::cout << externalPrecisions.format(CleanFmt) << std::endl;
			std::cout << "  Derivatives:" << std::endl;
			std::cout << externalDerivatives.format(CleanFmt) << std::endl;
		}
		if (externalPoint) {
			std::cout << " External seed:" << std::endl;
			std::cout << externalSeed.format(CleanFmt) << std::endl;
		}
		if (fitOK) {
			std::cout << " Fit results" << std::endl;
			std::cout << "  Parameters:" << std::endl;
			theVector.print();
			std::cout << "  Covariance matrix (bordered band part):"
					<< std::endl;
			theMatrix.printMatrix();
		}
	}
}

/// Print \link GblPoint GblPoints \endlink on trajectory
/**
 * \param [in] level print level (0: minimum, >0: more)
 */
void GblTrajectory::printPoints(unsigned int level) const {
	std::cout << "GblPoints " << std::endl;
	for (unsigned int iTraj = 0; iTraj < numTrajectories; ++iTraj) {
		std::vector<GblPoint>::const_iterator itPoint;
		for (itPoint = thePoints[iTraj].begin();
				itPoint < thePoints[iTraj].end(); ++itPoint) {
			itPoint->printPoint(level);
		}
	}
}

/// Print GblData blocks for trajectory
void GblTrajectory::printData() const {
	std::cout << "GblData blocks " << std::endl;
	std::vector<GblData>::const_iterator itData;
	for (itData = theData.begin(); itData < theData.end(); ++itData) {
		itData->printData();
	}
}

/// Get condition from band (decomposition)
/**
 * Return condition or 0. if band part is not positive definite
 * or -1. if the fit failed or has not been done yet.
 */
double GblTrajectory::getBandCondition() const {
	return (fitOK) ? theMatrix.getBandCondition() : -1.;
}

}
