/* Copyright 2011, Technische Universitaet Muenchen,
   Author: Karl Bicker

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
/** @addtogroup genfit */
/* @{ */

#ifndef GFTOOLS_H
#define GFTOOLS_H

#include <cmath>
#include <GFAbsTrackRep.h>
#include <GFDetPlane.h>
#include <GFException.h>
#include <GFTrack.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TDecompSVD.h>
//#include <TDecompBK.h>
//#include <TDecompChol.h>

/*! \namespace GFTools
 *  \brief Namespace for various tools, mainly smoothing.
 *
 * The GFTools namespace contains various functions, most of which er intended
 * to be used with the smoothing feature of the Kalman filter. These functions
 * allow to extract the smoothed information for any track which was fitted 
 * with smoothing enabled.
 */
namespace GFTools {
	
	/** @brief Get smoothed track position in plane coordinates
	 */
	TMatrixT<double> getSmoothedPos(GFTrack* trk, unsigned int irep, unsigned int ihit);

	/** @brief Get smoothed track covariance in plane coordinates
	 */
	TMatrixT<double> getSmoothedCov(GFTrack* trk, unsigned int irep, unsigned int ihit);

	/** @brief Get smoothed state vector and state covariance.
	 */
	bool getSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov);

	/** @brief Get smoothed state vector, state covariance and smoothing plane.
	 *
	 * The smoothed data returned by this function includes the detector plane in 
	 * which the data is calculated.
	 */
	bool getSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane);

	/** @brief Get smoothed state vector, state covariance and smoothing plane.
	 *
	 * The smoothed data returned by this function includes the detector plane in 
	 * which the data is calculated as well as the auxillary information for this
	 * plane.
	 * The smoothed covariance matrix at hit i is:
	 * \f[
	 * C_{k,\mbox{smooth}} = (C_{k,k-1,\rightarrow}^{-1}+C_{k,k+1,\leftarrow}^{-1})^{-1}. \\
	 * \f]
	 * The smoothed state vector at hit i is:
	 * \f[
	 * x_{k,\mbox{smooth}} = C_{k,\mbox{smooth}} \cdot ((C_{k,k-1,\rightarrow}^{-1}\cdot x_{k,k-1,\rightarrow} + C_{k,k+1,\leftarrow}^{-1}\cdot x_{k,k+1,\leftarrow})
	 * \f]
	 * The index k,k-1 denotes that the state vector or covariance matrix contains the information of all hits 1,2,....,k-1 and is extrapolated to hit k. The arrow pointing to the right indicates that the information is saved during the forward fit. The index k,k+1 indicates the the state vector or covariance matrix contains the information of all hits N,N-1,...,k+1 and is extrapolated to hit k. The left-poiting arrow indicates that the information is saved during the backward fit.
	 */
	bool getSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixT<double>& auxInfo);

	/** @brief Get smoothing plane.
	 */
	GFDetPlane getSmoothingPlane(GFTrack* trk, unsigned int irep, unsigned int ihit);

	/** @brief Get biased smoothed state vector and state covariance.
	 */
	bool getBiasedSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov);

	/** @brief Get biased smoothed state vector, state covariance and smoothing plane.
	 *
	 * The smoothed data returned by this function includes the detector plane in 
	 * which the data is calculated.
	 */
	bool getBiasedSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane);

	/** @brief Get biased smoothed state vector, state covariance and smoothing plane.
	 *
	 * The smoothed data returned by this function includes the detector plane in 
	 * which the data is calculated as well as the auxillary information for this
	 * plane.
	 */
	bool getBiasedSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixT<double>& auxInfo);

	/** @brief Get biased smoothed track position in plane coordinates
	 */
	TMatrixT<double> getBiasedSmoothedPos(GFTrack* trk, unsigned int irep, unsigned int ihit);

	/** @brief Get biased smoothed track covariance in plane coordinates
	 */
	TMatrixT<double> getBiasedSmoothedCov(GFTrack* trk, unsigned int irep, unsigned int ihit);


	/** @brief Invert a matrix, throwing GFException when inversion fails.
	 */
	void invertMatrix(const TMatrixT<double>& mat, TMatrixT<double>& inv);
	/** @brief Get smoothed chi2 for a specific hit (ihit).
	 *
	 * This function calculates a smoothed chi2 value for a specific hit from the
	 * the so called "biased" smoothed state vector and covariance matrix.
	 * If many values from different tracks in the same layer are calculated
	 * they should be chi2 distributed with dim(m) degree of freedom, where
	 * m is the measurement vector of the hit. So for a pixel detector that
	 * measures x and y dim(m) will be 2.
	 */
	double getSmoothedChiSqu(GFTrack* const trk, unsigned int irep, unsigned int ihit);

}

#endif

/** @} */ 
