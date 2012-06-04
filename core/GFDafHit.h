/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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
/** @addtogroup genfit
 * @{
 */

#ifndef GFDAFHIT_H
#define GFDAFHIT_H

#include<cmath>
#include<GFAbsRecoHit.h>
#include<GFAbsTrackRep.h>
#include<GFDetPlane.h>
#include<GFTools.h>
#include<stdlib.h>
#include<vector>

/** @brief Wrapper class for use with GFDaf.
 *
 * The GFDafHit is a hit class which acts as an effective hit. As the DAF is
 * capable of handling multiple hits in one plane, but GFKalman cannot handle
 * this, the GFDafHit combines all hits in one plane.
 */
class GFDafHit: public GFAbsRecoHit {
	public:

		GFDafHit(){};	

		~GFDafHit(){ fRawHits.clear(); };

		/** @brief Constructor adding hits.
		 *
		 * With this constructor, the GFDafHit is normally initialized. All the
		 * hits of the argument should be in one plane.
		 */
		GFDafHit(std::vector<GFAbsRecoHit*> HitsInPlane);

		/** @brief Get the measurement m,V
		 *
		 * Due to the nature of the GFDafHit, the coordinates returned are not 
		 * necessarily the coordinates of the real hits. There are two cases:
		 * if the GFDafHit contains only one hit, then the hit coordinates of
		 * this hit are returned. If however, the GFDafHit contains several hits,
		 * then the following formula is used to obtain the hit coordinates:
		 * \f[
		 * \mathbf{m} = \mathbf{V} \cdot \sum_{i}p_{i}\cdot \mathbf{V}_{i}^{-1}\cdot \mathbf{m}_{i}
		 * \f]
		 * with \f$\mathbf{m_{i}}\f$ the real hit coordinates, \f$\mathbf{V}_{i}\f$ 
		 * the real hit covariance and \f$p_{i}\f$ the weights. The sum runs over
		 * all real hits in the GFDafHit. \f$\mathbf{V}\f$ is the effective hit
		 * covariance as returned by getHitCov(). This calculation is only done,
		 * if the plane has changed from the last time gitHitCoord() was called.
		 *
		 * And now for the covariance matrix:
		 * There are two cases as well: if the GFDafHit contains
		 * only one hit, the covariance is calculated like:
		 * \f[
		 * \mathbf{V} = \frac{1}{p_{1}}\cdot\mathbf{V}_{1}
		 * \f]
		 * with the symbols analog to getHitCoord().
		 * If there are several hits in the GFDafHit, the following formula is
		 * used:
		 * \f[
		 * \mathbf{V} = \left( \sum_{i} p_{i} \cdot \left( \mathbf{V}_{i}\right)^{-1} \right)^{-1}
		 * \f]
		 * As before, these calculations are only done if the plane is different
		 * from the one getHitCov was last called.
		 */
		void getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TMatrixT<double>& statePred,const TMatrixT<double>& covPred,TMatrixT<double>& m, TMatrixT<double>& V);

		/** @brief Get the H matrix.
		 *
		 * Returns the H matrix of the first hit in the GFDafHit. This is valid
		 * because all hits are on the same plane.
		 */
		TMatrixT<double> getHMatrix(const GFAbsTrackRep* rep);

		GFDafHit* clone();

		/** @brief Get the detector plane.
		 *
		 * Returns the plane as returned by GFDetPlan() of the first hit in the
		 * GFDafHit. This is valid because all hits are on the same plane.
		 */
		const GFDetPlane& getDetPlane(GFAbsTrackRep* rep);

		/** @brief Set the weights.
		 *
		 * Set the weights as used by getHitCoord() and getHitCov().
		 */
		void setWeights(std::vector<double> weights);

		/** @brief Get the number of hits in the GFDafHit.
		 */
		unsigned int getNumHits() { return fRawHits.size(); };

		/** @brief Get at hit from the GFDafHit.
		 */
		GFAbsRecoHit* getHit(unsigned int ihit);

		/** @brief Get the name of the hit policy
		 *
		 * Returns the name of the hit policy of the first hit in the GFDafHit. This
		 * is valid because if there are several hits, they all have to be planar hits.
		 */
		const std::string& getPolicyName();

	private:
		bool fHitUpd;
		GFDetPlane fPl;
		std::vector<TMatrixT<double> > fCovInvs;
		std::vector<GFAbsRecoHit*> fRawHits;
		std::vector<double> fWeights;

	public:
		ClassDef(GFDafHit,2)

};
#endif

/** @} */

