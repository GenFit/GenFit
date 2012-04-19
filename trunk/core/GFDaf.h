/* Copyright 2011, Technische Universitaet Muenchen,
Authors: Karl Bicker, Christian Hoeppner

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

#ifndef GFDAF_H
#define GFDAF_H

#include<assert.h>
#include<cmath>
#include<GFAbsRecoHit.h>
#include<GFDafHit.h>
#include<GFKalman.h>
#include<GFTrack.h>
#include<stdlib.h>
#include<vector>

/** @brief Determinstic Annealing Filter (DAF) implementation. 
 *
 * @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 * @author Karl Bicker (Technische Universit&auml;t M&uuml;nchen)
 *
 * The DAF is an iterative Kalman filter with annealing. It is capable of 
 * fitting tracks which are contaminated with noise hits. The algorithm is 
 * taken from the references R. Fruehwirth & A. Strandlie, Computer Physics 
 * Communications 120 (1999) 197-214 and CERN thesis: Dissertation by Matthias 
 * Winkler.
 *
 * The weights which were assigned to the hits by the DAF are accessible by using the
 * bookkeeping object of the fitted track. The weight is stored under the key "dafWeight".
 * So to retrieve for example the weight of hit 10, fitted with track representation 2,
 * use GFTrack::getBK(2)->getNumber("dafWeight", 10,  double& wght).
 */
class GFDaf: GFKalman {
	public:

		GFDaf();
		~GFDaf() { };

		/** @brief Process a track using the DAF.
		 */
		void processTrack(GFTrack* trk);

		/** @brief Return the weights present after the track was processed.
		 *
		 * WARNING: This function is deprecated! Use the bookkeeping instead.
		 *		
		 * The DAF uses special effective hits defined in the class GFDafHit. A
		 * GFDafHit is a wrapper class and contains all the real hits from one plane.
		 * The structure of the return vector of getWeights allows to reconstruct in
		 * what way the hits were grouped: the outermost vector represents the track
		 * representation, there is one entry per track representation. The middle
		 * vector represents the effective hits, and the innermost vector contains
		 * the real hits contained in the corresponding effective hit.
		 */
		const std::vector<std::vector<std::vector<double> > > getWeights() { 
			std::cout<<"Warning: Using deprecated GFDaf::getWeights()! The weights of the hits are accessible in the bookkeeping of the track which was fitted, the key is \"dafWeight\""<<std::endl;
			return fWeights; 
		};

		/** @brief Set the probabilty cut for the weight calculation for the hits. 
		 *
		 * Currently supported are the values 0.01 0.005, and 0.001. The 
		 * corresponding chi2 cuts for different hits dimensionalities are hardcoded 
		 * in the implementation because I did not yet figure out how to calculate 
		 * them. Please feel very welcome to change the implementtion if you know how
		 * to do it.
		 */
		void setProbCut(double prob_cut);

		/** @brief Configure the annealing scheme.
		 *
		 * In the current implementation you need to provide at least one temperatures
		 * and not more then ten tempertatures.
		 */
		void setBetas(double b1,double b2=-1,double b3=-1.,double b4=-1.,double b5=-1.,double b6=-1.,double b7=-1.,double b8=-1.,double b9=-1.,double b10=-1.);

	private:

		/** @brief Initialize the GFDafHit and their weights before the fit.
		 */
		std::vector<GFDafHit*> initHitsWeights(GFTrack* trk);

		/** @brief Calculate the weights for the next fitting pass.
		  */
		std::vector<std::vector<double> > calcWeights(GFTrack* trk, double beta);

		/** @brief Copy the smoothing matrices from the source track to the target.
		 */
		void copySmoothing(GFTrack* source, GFTrack* target, int target_ire);

		void saveWeights(GFTrack* trk, const std::vector<std::vector<std::vector<double> > >& weights) const;

		std::vector<std::vector<std::vector<double> > > fWeights;
		std::vector<double> fBeta;
		std::map<int,double>  fchi2Cuts;

};

#endif

/** @} */

