/*
 * GblData.h
 *
 *  Created on: Aug 18, 2011
 *      Author: kleinwrt
 */

/** \file
 *  GblData definition.
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

#ifndef GBLDATA_H_
#define GBLDATA_H_

#include<iostream>
#include<vector>
#include <array>
#include<math.h>
#include "VMatrix.h"

#include "Eigen/Core"

//! Namespace for the general broken lines package
namespace gbl {

typedef Eigen::Matrix<double, 5, 5> Matrix5d;
typedef Eigen::Matrix<double, 4, 9> Matrix49d;

enum dataBlockType {
	None, InternalMeasurement, InternalKink, ExternalSeed, ExternalMeasurement
};

/// Data (block) for independent scalar measurement
/**
 * Data (block) containing value, precision and derivatives for measurements, kinks and seeds.
 * Created from attributes of GblPoints, used to construct linear equation system for track fit.
 */
class GblData {
public:
	GblData() : theLabel(0), theValue(0.), thePrecision(-1.), theDownWeight(0.), thePrediction(0.) {}; // not part of original GBL
	GblData(unsigned int aLabel, dataBlockType aType, double aValue,
			double aPrec, unsigned int aTraj = 0, unsigned int aPoint = 0,
			unsigned int aMeas = 0);
	GblData(const GblData&) = default;
	GblData& operator=(const GblData&) = default;
	GblData(GblData&&) = default;
	GblData& operator=(GblData&&) = default;
	virtual ~GblData();

	/// Add derivatives from measurement.
	/**
	 * Add (non-zero) derivatives to data block. Fill list of labels of used fit parameters.
	 * \tparam LocalDerivative  Local derivatives matrix
	 * \tparam ExtDerivative  External derivatives matrix
	 * \param [in] iRow Row index (0-4) in up to 5D measurement
	 * \param [in] labDer Labels for derivatives
	 * \param [in] matDer Derivatives (matrix) 'measurement vs track fit parameters'
	 * \param [in] iOff Offset for row index for additional parameters
	 * \param [in] derLocal Derivatives (matrix) for additional local parameters
	 * \param [in] extOff Offset for external parameters
	 * \param [in] extDer Derivatives for external Parameters
	 */
	template<typename LocalDerivative, typename ExtDerivative>
	void addDerivatives(unsigned int iRow,
			const std::array<unsigned int, 5> &labDer, const Matrix5d &matDer,
			unsigned int iOff,
			const Eigen::MatrixBase<LocalDerivative> &derLocal,
			unsigned int extOff,
			const Eigen::MatrixBase<ExtDerivative> &extDer);

	/// Add derivatives from kink.
	/**
	 * Add (non-zero) derivatives to data block. Fill list of labels of used fit parameters.
	 * \tparam ExtDerivative  External derivatives matrix
	 * \param [in] iRow Row index (0-3) in 2D kinks and steps
	 * \param [in] nDer Number of derivatives (7: kinks, 9: kinks+steps)
	 * \param [in] labDer Labels for derivatives
	 * \param [in] matDer Derivatives (matrix) 'kinks+steps vs track fit parameters'
	 * \param [in] extOff Offset for external parameters
	 * \param [in] extDer Derivatives for external Parameters
	 */
	template<typename ExtDerivative>
	void addDerivatives(unsigned int iRow, unsigned int nDer,
			const std::array<unsigned int, 9> &labDer, const Matrix49d &matDer,
			unsigned int extOff,
			const Eigen::MatrixBase<ExtDerivative> &extDer);

	void addDerivatives(const std::vector<unsigned int> &index,
			const std::vector<double> &derivatives);

	void setPrediction(const VVector &aVector);
	double setDownWeighting(unsigned int aMethod);
	double getChi2() const;
	void printData() const;
	unsigned int getLabel() const;
	dataBlockType getType() const;
	void getLocalData(double &aValue, double &aWeight, unsigned int &numLocal,
			unsigned int *&indLocal, double *&derLocal);
	void getAllData(double &aValue, double &aErr, unsigned int &numLocal,
			unsigned int *&indLocal, double *&derLocal, unsigned int &aTraj,
			unsigned int &aPoint, unsigned int &aMeas, unsigned int &aRow);
	void getResidual(double &aResidual, double &aVariance, double &aDownWeight,
			unsigned int &numLocal, unsigned int *&indLocal, double *&derLocal);
	void getResidual(double &aResidual, double &aVariance);

private:
	unsigned int theLabel; ///< Label (of corresponding point)
	unsigned int theRow; ///< Row number (of measurement)
	dataBlockType theType; ///< Type (None, InternalMeasurement, InternalKink, ExternalSeed, ExternalMeasurement)
	double theValue; ///< Value (residual)
	double thePrecision; ///< Precision (1/sigma**2)
	unsigned int theTrajectory; ///< Trajectory number
	unsigned int thePoint; ///< Point number (on trajectory)
	unsigned int theMeas; ///< Measurement number (at point)
	unsigned int theDWMethod; ///< Down-weighting method (0: None, 1: Tukey, 2: Huber, 3: Cauchy)
	double theDownWeight; ///< Down-weighting factor (0-1)
	double thePrediction; ///< Prediction from fit
	// standard local parameters (curvature, offsets), fixed size
	unsigned int theNumLocal; ///< Number of (non zero) local derivatives (max 9 for kinks+steps)
	unsigned int theParameters[9]; ///< List of parameters (with non zero derivatives)
	double theDerivatives[9]; ///< List of derivatives for fit
	// more local parameters, dynamic size
	std::vector<unsigned int> moreParameters; ///< List of fit parameters (with non zero derivatives)
	std::vector<double> moreDerivatives; ///< List of derivatives for fit
};

template<typename LocalDerivative, typename ExtDerivative>
void GblData::addDerivatives(unsigned int iRow,
		const std::array<unsigned int, 5> &labDer, const Matrix5d &matDer,
		unsigned int iOff, const Eigen::MatrixBase<LocalDerivative> &derLocal,
		unsigned int extOff, const Eigen::MatrixBase<ExtDerivative> &extDer) {

	unsigned int nParMax = 5 + derLocal.cols() + extDer.cols();
	theRow = iRow - iOff;
	if (nParMax > 9) {
		// dynamic data block size
		moreParameters.reserve(nParMax); // have to be sorted
		moreDerivatives.reserve(nParMax);

		for (int i = 0; i < derLocal.cols(); ++i) // local derivatives
				{
			if (derLocal(iRow - iOff, i)) {
				moreParameters.push_back(i + 1);
				moreDerivatives.push_back(derLocal(iRow - iOff, i));
			}
		}

		for (int i = 0; i < extDer.cols(); ++i) // external derivatives
				{
			if (extDer(iRow - iOff, i)) {
				moreParameters.push_back(extOff + i + 1);
				moreDerivatives.push_back(extDer(iRow - iOff, i));
			}
		}

		for (size_t i = 0; i < labDer.size(); ++i) // curvature, offset derivatives
				{
			if (labDer[i] and matDer(iRow, i)) {
				moreParameters.push_back(labDer[i]);
				moreDerivatives.push_back(matDer(iRow, i));
			}
		}
	} else {
		// simple (static)  data block
		for (int i = 0; i < derLocal.cols(); ++i) // local derivatives
				{
			if (derLocal(iRow - iOff, i)) {
				theParameters[theNumLocal] = i + 1;
				theDerivatives[theNumLocal] = derLocal(iRow - iOff, i);
				theNumLocal++;
			}
		}

		for (int i = 0; i < extDer.cols(); ++i) // external derivatives
				{
			if (extDer(iRow - iOff, i)) {
				theParameters[theNumLocal] = extOff + i + 1;
				theDerivatives[theNumLocal] = extDer(iRow - iOff, i);
				theNumLocal++;
			}
		}
		for (size_t i = 0; i < labDer.size(); ++i) // curvature, offset derivatives
				{
			if (labDer[i] and matDer(iRow, i)) {
				theParameters[theNumLocal] = labDer[i];
				theDerivatives[theNumLocal] = matDer(iRow, i);
				theNumLocal++;
			}
		}
	}
}

template<typename ExtDerivative>
void GblData::addDerivatives(unsigned int iRow,unsigned int nDer,
		const std::array<unsigned int, 9> &labDer, const Matrix49d &matDer,
		unsigned int extOff, const Eigen::MatrixBase<ExtDerivative> &extDer) {

	unsigned int nParMax = nDer + extDer.cols();
	theRow = iRow;
	if (nParMax > 9) {
		// dynamic data block size
		moreParameters.reserve(nParMax); // have to be sorted
		moreDerivatives.reserve(nParMax);

		for (int i = 0; i < extDer.cols(); ++i) // external derivatives
				{
			if (extDer(iRow, i)) {
				moreParameters.push_back(extOff + i + 1);
				moreDerivatives.push_back(extDer(iRow, i));
			}
		}

		for (size_t i = 0; i < nDer; ++i) // curvature, offset derivatives
				{
			if (labDer[i] and matDer(iRow, i)) {
				moreParameters.push_back(labDer[i]);
				moreDerivatives.push_back(matDer(iRow, i));
			}
		}
	} else {
		// simple (static) data block
		for (size_t i = 0; i < nDer; ++i) // curvature, offset derivatives
				{
			if (labDer[i] and matDer(iRow, i)) {
				theParameters[theNumLocal] = labDer[i];
				theDerivatives[theNumLocal] = matDer(iRow, i);
				theNumLocal++;
			}
		}
	}
}

}

#endif /* GBLDATA_H_ */
