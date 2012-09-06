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

#ifndef GFKALMAN_H
#define GFKALMAN_H

#include <map>
#include <iostream>

#include "TMatrixT.h"


class GFAbsRecoHit;
class GFAbsTrackRep;
class GFTrack;

/** @brief Generic Kalman Filter implementation
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * The Kalman Filter operates on genfit GFTrack objects. It is an implementation
 * of the Kalman Filter algebra that uses the genfit interface classes 
 * GFAbsRecoHit and GFAbsTrackRep in order to be independent from the specific
 * detector geometry and the details of the track parameterization /
 * track extraoplation engine.
 *
 * The Kalman Filter can use hits from several detectors in a single fit 
 * to estimate the parameters of several track representations in parallel.
 */
class GFKalman {
public:

  //friend class KalmanTester; // gives the Tester access to private methods

  // Constructors/Destructors ---------
  GFKalman();
  ~GFKalman();

  // Operators
  /** @brief Operator for use with STL.
   *
   * This operator allows to use the std::foreach algorithm with an
   * STL container o GFTrack* objects.
   */
  inline void operator()(GFTrack* track){processTrack(track);}
  
  /** @brief Operator for use with STL.
   *
   * This operator allows to use the std::foreach algorithm with an
   * STL container o GFTrack* objects.
   */
  inline void operator()(std::pair<int,GFTrack*> tr){processTrack(tr.second);}

  // Operations ----------------------

  /** @brief Switch lazy error handling.
   *
   * This is a historically left-over method and shall be deleted some time
   */
  void setLazy(Int_t flag){
    static_cast<void>(flag);
    std::cerr<<"Using outdates setLazy method of class GFKalman:"<<std::endl;
  }

  /** @brief Set number of iterations for Kalman Filter
   *
   * One iteration is one forward pass plus one backward pass
   */
  void setNumIterations(Int_t i){fNumIt=i;}

  /** @brief Performs fit on a GFTrack.
   *
   * The hits are processed in the order in which they are stored in the GFTrack
   * object. Sorting of hits in space has to be done before!
   */
  void processTrack(GFTrack* trk);

  /** @brief Performs fit on a GFTrack beginning with the current hit.
   */
  void fittingPass(GFTrack*,int dir); // continues track from lastHitInFit

  /** @brief Calculates chi2 of a given hit with respect to a 
   * given track representation.
   */
  double getChi2Hit(GFAbsRecoHit*, GFAbsTrackRep*);

  /** @brief Sets the inital direction of the track fit (1 for inner to outer,
   * or -1 for outer to inner). The standard is 1 and is set in the ctor
   */
  void setInitialDirection(int d){fInitialDirection=d;}

  /** @brief Set the blowup factor (see blowUpCovs() )
   */
  void setBlowUpFactor(double f){fBlowUpFactor=f;}
  
  // Protected Methods -----------------
protected:  
  /** @brief this is needed to blow up the covariance matrix before a fitting pass
   * drops off-diagonal elements and blows up diagonal by blowUpFactor
   */
  void blowUpCovs(GFTrack* trk);

  // Private Methods -----------------
private:
  /** @brief One Kalman step.
   *
   * Performs
   * - Extrapolation to detector plane of the hit
   * - Calculation of residual and Kalman Gain
   * - Update of track representation state and chi2
   *
   */
  void processHit(GFTrack*, int, int, int);
  
  /** @brief Used to switch between forward and backward filtering
   */
  void switchDirection(GFTrack* trk); // switches the direction of propagation for all reps

  /** @brief Calculate Kalman Gain
   */
  TMatrixT<double> calcGain(const TMatrixT<double>& cov, 
						const TMatrixT<double>& HitCov,
						const TMatrixT<double>& H);

  /** @brief this returns the reduced chi2 increment for a hit
   */
  double chi2Increment(const TMatrixT<double>& r,const TMatrixT<double>& H,
		       const TMatrixT<double>& cov,const TMatrixT<double>& V);



  int fInitialDirection;
  Int_t fNumIt;
  double fBlowUpFactor;
  bool fSmooth;
  bool fSmoothFast;




};


#endif

/** @} */
