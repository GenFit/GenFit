/*
 * GblUtilities.h
 *
 *  Created on: 6 Nov 2018
 *      Author: kleinwrt
 */

/** \file
 *  Definitions for GBL utilities.
 *
 *  \author Claus Kleinwort, DESY, 2018 (Claus.Kleinwort@desy.de)
 *
 *  \copyright
 *  Copyright (c) 2018-2024 Deutsches Elektronen-Synchroton,
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

#ifndef GBLUTILITIES_H_
#define GBLUTILITIES_H_

#include "Eigen/Dense"
#include<iostream>

//! Namespace for the general broken lines package
namespace gbl {

double unrm();
double unif();
Eigen::Matrix<double, 5, 5> gblSimpleJacobian(double ds, double cosl,
		double bfac);
double gblMultipleScatteringError(double qbyp, double xbyx0);

/// Prediction on helix
/**
 * Prediction at intersection of helix and measurement plane.
 */
class GblHelixPrediction {
public:
	GblHelixPrediction(double sArc, const Eigen::Vector2d &aPred,
			const Eigen::Vector3d &tDir, const Eigen::Vector3d &uDir,
			const Eigen::Vector3d &vDir, const Eigen::Vector3d &nDir,
			const Eigen::Vector3d &aPos);
	virtual ~GblHelixPrediction();
	double getArcLength() const;
	const Eigen::Vector2d& getMeasPred() const;
	const Eigen::Vector3d& getPosition() const;
	const Eigen::Vector3d& getDirection() const;
	double getCosIncidence() const;
	Eigen::Matrix<double, 2, 3> getCurvilinearDirs() const;

private:
	const double sarc; ///< arc-length at prediction
	const Eigen::Vector2d pred; ///< prediction for measurement (u,v)
	const Eigen::Vector3d tdir; ///< track direction at prediction
	const Eigen::Vector3d udir; ///< measurement direction for u
	const Eigen::Vector3d vdir; ///< measurement direction for v
	const Eigen::Vector3d ndir; ///< normal to measurement plane
	const Eigen::Vector3d pos; ///< position at prediction
	Eigen::Matrix3d global2meas; ///< transformation into measurement system
};

///Simple helix
/**
 * Circle in XY plane, straight line in ZS.
 */
class GblSimpleHelix {
public:
	GblSimpleHelix(double aRinv, double aPhi0, double aDca, double aDzds,
			double aZ0);
	virtual ~GblSimpleHelix();
	double getPhi(double aRadius) const;
	double getArcLengthR(double aRadius) const;
	double getArcLengthXY(double xPos, double yPos) const;
	void moveToXY(double xPos, double yPos, double &newPhi0, double &newDca,
			double &newZ0) const;
	GblHelixPrediction getPrediction(const Eigen::Vector3d &refPos,
			const Eigen::Vector3d &uDir, const Eigen::Vector3d &vDir) const;

private:
	const double rinv; ///< curvature (1/Radius)
	const double phi0; ///< azimuth at PCA (point of closest approach to origin in XY plane, defines arc-length S=0)
	const double dca; ///< distance to origin in XY plane at PCA
	const double dzds; ///< slope in ZS plane (dZ/dS)
	const double z0; ///< offset in ZS plane
	const double cosPhi0; ///< cos(phi0)
	const double sinPhi0; ///< sin(phi0)
	const double xRelCenter; ///< X position of circle center / R
	const double yRelCenter; ///< Y position of circle center / R
};

/// Detector layer
/**
 * Alignable (rigid body) planar detector layer.
 */
class GblDetectorLayer {
public:
	GblDetectorLayer(const std::string aName, const unsigned int aLayer,
			const int aDim, const double thickness, Eigen::Vector3d &aCenter,
			Eigen::Vector2d &aResolution, Eigen::Vector2d &aPrecision,
			Eigen::Matrix3d &measTrafo, Eigen::Matrix3d &alignTrafo);
	virtual ~GblDetectorLayer();
	void print() const;
	void printMP2Constraint() const;
	unsigned int getRigidBodyGlobalLabel(const unsigned int aPar) const;
	unsigned int getLayerID() const;
	double getRadiationLength() const;
	Eigen::Vector2d getResolution() const;
	Eigen::Vector2d getPrecision() const;
	Eigen::Vector3d getCenter() const;
	Eigen::Matrix3d getMeasSystemDirs() const;
	Eigen::Matrix3d getAlignSystemDirs() const;
	GblHelixPrediction intersectWithHelix(GblSimpleHelix hlx) const;
	Eigen::Matrix<double, 3, 6> getRigidBodyDerGlobal(Eigen::Vector3d &position,
			Eigen::Vector3d &direction) const;
	Eigen::Matrix<double, 2, 6> getRigidBodyDerLocal(Eigen::Vector3d &position,
			Eigen::Vector3d &direction) const;
	Eigen::Matrix<double, 6, 6> getTrafoGlobalToLocal(Eigen::Vector3d &offset,
			Eigen::Matrix3d &rotation) const;
	Eigen::Matrix<double, 6, 6> getTrafoLocalToGlobal(Eigen::Vector3d &offset,
			Eigen::Matrix3d &rotation) const;

private:
	std::string name; ///< name
	unsigned int layer; ///< layer ID
	unsigned int measDim; ///< measurement dimension (1 or 2)
	double xbyx0; ///< normalized material thickness
	Eigen::Vector3d center; ///< center
	Eigen::Vector2d resolution; ///< measurements resolution
	Eigen::Vector2d precision; ///< measurements precision
	Eigen::Vector3d udir; ///< 1. measurement direction
	Eigen::Vector3d vdir; ///< 2. measurement direction
	Eigen::Vector3d ndir; ///< normal to measurement plane
	Eigen::Matrix3d global2meas; ///< transformation into measurement system
	Eigen::Matrix3d global2align; ///< transformation into (local) alignment system
	bool alignInMeasSys; ///< alignment == measurement system?
};

}
#endif /* GBLUTILITIES_H_ */
