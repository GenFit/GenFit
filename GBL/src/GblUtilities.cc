/*
 * exampleUtil.cpp
 *
 *  Created on: 6 Nov 2018
 *      Author: kleinwrt
 */

/** \file
 *  Utilities for GBL applications.
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

#include "GblUtilities.h"

using namespace Eigen;

namespace gbl {

typedef Eigen::Matrix<double, 5, 5> Matrix5d;

/// Multiple scattering error
/**
 * Angular error in plane, simple model (Rossi, Greisen, (1941))
 * \param [in] qbyp    q/p [1/GeV]
 * \param [in] xbyx0   thickness / radiation length
 */
double gblMultipleScatteringError(double qbyp, double xbyx0) {
	return 0.015 * fabs(qbyp) * sqrt(xbyx0);
}

/// Simple jacobian
/**
 * Simple jacobian for (q/p, slopes, offsets) in curvilinear system,
 * constant magnetic field in Z direction, quadratic in arc length difference.
 *
 * \param [in] ds    arc-length
 * \param [in] cosl  cos(lambda)
 * \param [in] bfac  Bz*c
 * \return jacobian to move by 'ds' on trajectory
 */
Matrix5d gblSimpleJacobian(double ds, double cosl, double bfac) {
	Matrix5d jac;
	jac.setIdentity();
	jac(1, 0) = -bfac * ds * cosl;
	jac(3, 0) = -0.5 * bfac * ds * ds * cosl;
	jac(3, 1) = ds;
	jac(4, 2) = ds;
	return jac;
}

///  unit normal distribution, Box-Muller method, polar form
double unrm() {
	static double unrm2 = 0.0;
	static bool cached = false;
	if (!cached) {
		double x, y, r;
		do {
			x = 2.0 * static_cast<double>(rand())
					/ static_cast<double>(RAND_MAX) - 1;
			y = 2.0 * static_cast<double>(rand())
					/ static_cast<double>(RAND_MAX) - 1;
			r = x * x + y * y;
		} while (r == 0.0 || r > 1.0);
		// (x,y) in unit circle
		double d = sqrt(-2.0 * log(r) / r);
		double unrm1 = x * d;
		unrm2 = y * d;
		cached = true;
		return unrm1;
	} else {
		cached = false;
		return unrm2;
	}
}

///  uniform distribution [0..1]
double unif() {
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

/// Create helix prediction.
/**
 * Prediction at intersection of helix and measurement plane.
 *
 * \param [in] sArc     arc length
 * \param [in] aPred    prediction for measurement (u,v)
 * \param [in] tDir     track direction at prediction
 * \param [in] uDir     measurement direction for u
 * \param [in] vDir     measurement direction for v
 * \param [in] nDir     normal to measurement plane
 * \param [in] aPos     position at prediction
 */
GblHelixPrediction::GblHelixPrediction(double sArc, const Vector2d &aPred,
		const Vector3d &tDir, const Vector3d &uDir, const Vector3d &vDir,
		const Vector3d &nDir, const Vector3d &aPos) :
		sarc(sArc), pred(aPred), tdir(tDir), udir(uDir), vdir(vDir), ndir(nDir), pos(
				aPos) {
	global2meas << uDir.transpose(), vDir.transpose(), nDir.transpose();
}

GblHelixPrediction::~GblHelixPrediction() {
}

/// Get arc-length.
double GblHelixPrediction::getArcLength() const {
	return sarc;
}

/// Get (measurement) prediction.
const Vector2d& GblHelixPrediction::getMeasPred() const {
	return pred;
}

/// Get position.
const Vector3d& GblHelixPrediction::getPosition() const {
	return pos;
}

/// Get position.
const Vector3d& GblHelixPrediction::getDirection() const {
	return tdir;
}

/// Get cosine of incidence
double GblHelixPrediction::getCosIncidence() const {
	return tdir.dot(ndir);
}

/// Get curvilinear directions (U,V)
/*
 * Curvilinear system: track direction T, U = Z x T / |Z x T|, V = T x U
 */
Eigen::Matrix<double, 2, 3> GblHelixPrediction::getCurvilinearDirs() const {
	const double cosTheta = tdir[2];
	const double sinTheta = sqrt(tdir[0] * tdir[0] + tdir[1] * tdir[1]);
	const double cosPhi = tdir[0] / sinTheta;
	const double sinPhi = tdir[1] / sinTheta;
	Eigen::Matrix<double, 2, 3> uvDir;
	uvDir << -sinPhi, cosPhi, 0., -cosPhi * cosTheta, -sinPhi * cosTheta, sinTheta;
	return uvDir;
}

/// Create simple helix.
/**
 * Helix for constant magnetic field in Z direction.
 *
 * \param [in] aRinv      curvature (1/R)
 * \param [in] aPhi0      azimuth at PCA
 * \param [in] aDca       XY distance at PCA
 * \param [in] aDzds      slope in ZS (tanLambda)
 * \param [in] aZ0        offset in ZS
 */
GblSimpleHelix::GblSimpleHelix(double aRinv, double aPhi0, double aDca,
		double aDzds, double aZ0) :
		rinv(aRinv), phi0(aPhi0), dca(aDca), dzds(aDzds), z0(aZ0), cosPhi0(
				cos(phi0)), sinPhi0(sin(phi0)), xRelCenter(
				-(1. - dca * rinv) * sinPhi0), yRelCenter(
				(1. - dca * rinv) * cosPhi0) {
}

GblSimpleHelix::~GblSimpleHelix() {
}

/// Get phi (of point on circle) for given radius (to ref. point)
/**
 * ( |dca| < radius < |rad-2*dca|, from H1/cjfphi, not restricted to -Pi .. +Pi )
 *
 * \param[in] aRadius radius
 */
double GblSimpleHelix::getPhi(double aRadius) const {
	double arg = (0.5 * rinv * (aRadius * aRadius + dca * dca) - dca)
			/ (aRadius * (1.0 - rinv * dca));
	return asin(arg) + phi0;
}

/// Get (2D) arc length for given radius (to ref. point)
/**
 * ( |dca| < radius < |rad-2*dca|, from H1/cjfsxy )
 *
 * \param[in] aRadius radius
 */
double GblSimpleHelix::getArcLengthR(double aRadius) const {
	double arg = (0.5 * rinv * (aRadius * aRadius + dca * dca) - dca)
			/ (aRadius * (1.0 - rinv * dca));
	if (fabs(arg) >= 1.) {
		std::cout << " bad arc " << aRadius << " " << rinv << " " << dca
				<< std::endl;
		return 0.;
	}
	// line
	if (rinv == 0)
		return sqrt(aRadius * aRadius - dca * dca);
	// helix
	double sxy = asin(aRadius * rinv * sqrt(1.0 - arg * arg)) / rinv;
	if (0.5 * rinv * rinv * (aRadius * aRadius - dca * dca) - 1. + rinv * dca
			> 0.)
		sxy = M_PI / fabs(rinv) - sxy;
	return sxy;
}

/// Get (2D) arc length for given point.
/**
 * \param [in] xPos   X Position
 * \param [in] yPos   Y Position
 */
double GblSimpleHelix::getArcLengthXY(double xPos, double yPos) const {
// line
	if (rinv == 0)
		return cosPhi0 * xPos + sinPhi0 * yPos;
// helix
	double dx = xPos * rinv - xRelCenter;
	double dy = yPos * rinv - yRelCenter;
	double dphi = atan2(dx, -dy) - phi0;
	if (dphi > M_PI)
		dphi -= 2.0 * M_PI;
	else if (dphi < -M_PI)
		dphi += 2.0 * M_PI;
	return dphi / rinv;
}

/// Move to new reference point (X,Y)
/**
 * \param [in] xPos      X Position
 * \param [in] yPos      Y Position
 * \param [out] newPhi0  new phi0
 * \param [out] newDca   new dca
 * \param [out] newZ0    new z0
 */
void GblSimpleHelix::moveToXY(double xPos, double yPos, double &newPhi0,
		double &newDca, double &newZ0) const {
// start values
	newPhi0 = phi0;
	newDca = dca;
	newZ0 = z0;
// Based on V. Karimaki, NIM A305 (1991) 187-191, eqn (19)
	const double u = 1. - rinv * dca;
	const double dp = -xPos * sinPhi0 + yPos * cosPhi0 + dca;
	const double dl = xPos * cosPhi0 + yPos * sinPhi0;
	const double sa = 2. * dp - rinv * (dp * dp + dl * dl);
	const double sb = rinv * xPos + u * sinPhi0;
	const double sc = -rinv * yPos + u * cosPhi0;
	const double sd = sqrt(1. - rinv * sa);
// transformed parameters
	double sArc;
	if (rinv == 0.) {
		newDca = dp;
		sArc = dl;
	} else {
		newPhi0 = atan2(sb, sc);
		newDca = sa / (1. + sd);
		double dphi = newPhi0 - phi0;
		if (dphi > M_PI)
			dphi -= 2.0 * M_PI;
		else if (dphi < -M_PI)
			dphi += 2.0 * M_PI;
		sArc = dphi / rinv;
	}
	newZ0 += sArc * dzds;
}

/// Get prediction
/*
 * \param [in] refPos  reference position on detector plane
 * \param [in] uDir    measurement direction 'u'
 * \param [in] vDir    measurement direction 'v'
 */
GblHelixPrediction GblSimpleHelix::getPrediction(const Eigen::Vector3d &refPos,
		const Eigen::Vector3d &uDir, const Eigen::Vector3d &vDir) const {
// normal to (u,v) measurement plane
	Vector3d nDir = uDir.cross(vDir).normalized();
// ZS direction
	const double cosLambda = 1. / sqrt(1. + dzds * dzds);
	const double sinLambda = dzds * cosLambda;
	double sArc2D;
	Vector3d dist, pos, tDir;
// line (or helix)
	if (rinv == 0.) {
		// track direction
		tDir << cosLambda * cosPhi0, cosLambda * sinPhi0, sinLambda;
		// distance (of point at dca to reference)
		Vector3d pca(dca * sinPhi0, -dca * cosPhi0, z0);
		dist = pca - refPos;
		// arc-length
		double sArc3D = -dist.dot(nDir) / tDir.dot(nDir);
		sArc2D = sArc3D * cosLambda;
		// position at prediction
		pos = pca + sArc3D * tDir;
		// distance (of point at sArc to reference)
		dist = pos - refPos;
	} else {
		// initial guess of 2D arc-length
		sArc2D = this->getArcLengthXY(refPos(0), refPos(1));
		unsigned int nIter = 0;
		while (nIter < 10) {
			nIter += 1;
			// track direction
			const double dPhi = sArc2D * rinv;
			const double cosPhi = cos(phi0 + dPhi);
			const double sinPhi = sin(phi0 + dPhi);
			tDir << cosLambda * cosPhi, cosLambda * sinPhi, sinLambda;
			// position at prediction
			pos << (xRelCenter + sinPhi) / rinv, (yRelCenter - cosPhi) / rinv, z0
					+ dzds * sArc2D;
			// distance (of point at sArc to reference)
			dist = pos - refPos;
			// arc-length correction (linearizing helix at sArc)
			const double sCorr3D = -dist.dot(nDir) / tDir.dot(nDir);
			if (fabs(sCorr3D) > 0.00001) {
				// iterate
				sArc2D += sCorr3D * cosLambda;
			} else {
				// converged
				break;
			}
		}
	}
// projections on measurement directions
	Vector2d pred(dist.dot(uDir), dist.dot(vDir));
	return GblHelixPrediction(sArc2D, pred, tDir, uDir, vDir, nDir, pos);
}

/// Create a detector layer.
/**
 * Create planar detector layer with 1D or 2D measurement (u,v).
 *
 * \param [in] aName          name
 * \param [in] aLayer         layer ID
 * \param [in] aDim           dimension (1,2)
 * \param [in] thickness      thickness / radiation_length
 * \param [in] aCenter        center of detector (origin of local systems)
 * \param [in] aResolution    resolution vector
 * \param [in] aPrecision     diagonal of precision matrix
 * \param [in] measTrafo      matrix of row vectors defining local measurement system
 * \param [in] alignTrafo     matrix of row vectors defining local alignment system
 */
GblDetectorLayer::GblDetectorLayer(const std::string aName,
		const unsigned int aLayer, const int aDim, const double thickness,
		Eigen::Vector3d &aCenter, Eigen::Vector2d &aResolution,
		Eigen::Vector2d &aPrecision, Eigen::Matrix3d &measTrafo,
		Eigen::Matrix3d &alignTrafo) :
		name(aName), layer(aLayer), measDim(aDim), xbyx0(thickness), center(
				aCenter), resolution(aResolution), precision(aPrecision), global2meas(
				measTrafo), global2align(alignTrafo) {
	udir = global2meas.row(0);
	vdir = global2meas.row(1);
	ndir = global2meas.row(2);
	alignInMeasSys = global2meas.isApprox(global2align);
}

GblDetectorLayer::~GblDetectorLayer() {
}

/// Print GblDetectorLayer.
void GblDetectorLayer::print() const {
	IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::cout << " Layer " << name << " " << layer << " : " << measDim << "D, "
			<< xbyx0 << " X0, @ " << center.transpose().format(CleanFmt)
			<< ", res " << resolution.transpose().format(CleanFmt) << ", udir "
			<< udir.transpose().format(CleanFmt) << ", vdir "
			<< vdir.transpose().format(CleanFmt) << std::endl;
}

/// Print MP2 constraint.
/*
 * Alignment for **single** 1D measurement outside measurement system requires constraint (for offsets in v direction).
 * If there are multiple 1D measurements for an alignable ('layer') with different orientations the corresponding
 * constraints must be ignored.
 */
void GblDetectorLayer::printMP2Constraint() const {
	if (measDim > 1 or alignInMeasSys)
		return;
	// transform vdir into alignment system
	Eigen::Vector3d unMeasured = global2align * vdir;
	std::cout << "Constraint 0. ! fix unmeasured direction in " << name
			<< std::endl;
	for (int p = 0; p < 3; p++) {
		// 'zero' suppression
		if (fabs(unMeasured(p)) > 1.0e-10)
			std::cout << " " << getRigidBodyGlobalLabel(p) << " "
					<< unMeasured(p) << std::endl;
	}
}

/// Get global label
/*
 * Get global label for rigid body alignment parameters
 * (3 offsets, 3 rotations) in (local) alignment system.
 *
 * \param[in] aPar  parameter index (0..5)
 */
unsigned int GblDetectorLayer::getRigidBodyGlobalLabel(
		const unsigned int aPar) const {
	return layer * 10 + aPar + 1;
}

/// Get layer ID
unsigned int GblDetectorLayer::getLayerID() const {
	return layer;
}

/// Get radiation length.
double GblDetectorLayer::getRadiationLength() const {
	return xbyx0;
}

/// Get resolution.
Eigen::Vector2d GblDetectorLayer::getResolution() const {
	return resolution;
}

/// Get precision.
Eigen::Vector2d GblDetectorLayer::getPrecision() const {
	return precision;
}

/// Get center.
Eigen::Vector3d GblDetectorLayer::getCenter() const {
	return center;
}

/// Get directions of measurement system.
/**
 * Matrix from row vectors (transformation from global to measurement system)
 */
Eigen::Matrix3d GblDetectorLayer::getMeasSystemDirs() const {
	return global2meas;;
}

/// Get directions of alignment system.
/**
 * Matrix from row vectors (transformation from global to alignment system)
 */
Eigen::Matrix3d GblDetectorLayer::getAlignSystemDirs() const {
	return global2align;;
}

/// Intersect with helix.
/**
 * \param [in]  hlx  helix
 */
GblHelixPrediction GblDetectorLayer::intersectWithHelix(
		GblSimpleHelix hlx) const {
	return hlx.getPrediction(center, udir, vdir);
}

/// Get rigid body derivatives in global frame.
/**
 * \param[in] position   position (of prediction or measurement)
 * \param[in] direction  track direction
 */
Matrix<double, 3, 6> GblDetectorLayer::getRigidBodyDerGlobal(
		Eigen::Vector3d &position, Eigen::Vector3d &direction) const {
// lever arms (for rotations)
	Vector3d dist = position;
// dr/dm (residual vs measurement, 1-tdir*ndir^t/tdir*ndir)
	Matrix3d drdm = Matrix3d::Identity()
			- (direction * ndir.transpose()) / (direction.transpose() * ndir);
// dm/dg (measurement vs 6 rigid body parameters, global system)
	Matrix<double, 3, 6> dmdg = Matrix<double, 3, 6>::Zero();
	dmdg(0, 0) = 1.;
	dmdg(0, 4) = -dist(2);
	dmdg(0, 5) = dist(1);
	dmdg(1, 1) = 1.;
	dmdg(1, 3) = dist(2);
	dmdg(1, 5) = -dist(0);
	dmdg(2, 2) = 1.;
	dmdg(2, 3) = -dist(1);
	dmdg(2, 4) = dist(0);
// drl/dg (local residuals vs rigid body parameters)
	return global2meas * drdm * dmdg;
}

/// Get rigid body derivatives in local (alignment) frame (rotated in measurement plane).
/**
 * The orthogonal alignment frame differs from measurement frame only by rotations
 * around normal to measurement plane.
 *
 * Equivalent to:
 * \code
 * getRigidBodyDerGlobal(position, direction) * getTrafoLocalToGlobal(center, global2align)
 * \endcode
 *
 * \param[in] position   position (of prediction or measurement)
 * \param[in] direction  track direction
 */
Matrix<double, 2, 6> GblDetectorLayer::getRigidBodyDerLocal(
		Eigen::Vector3d &position, Eigen::Vector3d &direction) const {
	// track direction in local system
	Vector3d tLoc = global2align * direction;
	// local slopes
	const double uSlope = tLoc[0] / tLoc[2];
	const double vSlope = tLoc[1] / tLoc[2];
	// lever arms (for rotations)
	Vector3d dist = global2align * (position - center);
	const double uPos = dist[0];
	const double vPos = dist[1];
	// wPos = 0 (in detector plane)
	// drl/dg (local residuals (du,dv) vs rigid body parameters)
	Matrix<double, 2, 6> drldg;
	drldg << 1.0, 0.0, -uSlope, vPos * uSlope, -uPos * uSlope, vPos, 0.0, 1.0, -vSlope, vPos
			* vSlope, -uPos * vSlope, -uPos;
	// avoid numerics in case of unit transformation (below)
	if (alignInMeasSys)
		return drldg;
	// local (alignment) to measurement system
	Matrix3d local2meas = global2meas * global2align.transpose();
	return local2meas.block<2, 2>(0, 0) * drldg;
}

/// Get transformation for rigid body derivatives from global to local (alignment) system.
/**
 * local = rotation * (global-offset)
 *
 * \param[in] offset    offset of alignment system
 * \param[in] rotation  rotation of alignment system
 */
Matrix<double, 6, 6> GblDetectorLayer::getTrafoGlobalToLocal(
		Eigen::Vector3d &offset, Eigen::Matrix3d &rotation) const {
	// transformation global to local
	Matrix<double, 6, 6> glo2loc = Matrix<double, 6, 6>::Zero();
	Matrix3d leverArms;
	leverArms << 0., offset[2], -offset[1], -offset[2], 0., offset[0], offset[1], -offset[0], 0.;
	glo2loc.block<3, 3>(0, 0) = rotation;
	glo2loc.block<3, 3>(0, 3) = -rotation * leverArms;
	glo2loc.block<3, 3>(3, 3) = rotation;
	return glo2loc;
}

/// Get transformation for rigid body derivatives from local (alignment) to global system.
/**
 * local = rotation * (global-offset)
 *
 * \param[in] offset    offset of alignment system
 * \param[in] rotation  rotation of alignment system
 */
Matrix<double, 6, 6> GblDetectorLayer::getTrafoLocalToGlobal(
		Eigen::Vector3d &offset, Eigen::Matrix3d &rotation) const {
	// transformation local to global
	Matrix<double, 6, 6> loc2glo = Matrix<double, 6, 6>::Zero();
	Matrix3d leverArms;
	leverArms << 0., offset[2], -offset[1], -offset[2], 0., offset[0], offset[1], -offset[0], 0.;
	loc2glo.block<3, 3>(0, 0) = rotation.transpose();
	loc2glo.block<3, 3>(0, 3) = leverArms * rotation.transpose();
	loc2glo.block<3, 3>(3, 3) = rotation.transpose();
	return loc2glo;
}

}
