/* Copyright 2013, Ludwig-Maximilians Universität München,
   Authors: Tobias Schlüter & Johannes Rauch

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

#ifndef genfit_DAF_h
#define genfit_DAF_h

#include "AbsKalmanFitter.h"

#include <vector>
#include <map>
#include <memory>


namespace genfit {

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
 * The weights which were assigned to the hits by the DAF are accessible in the MeasurementOnPlane objects
 * in the KalmanFitterInfo objects.
 */
class DAF : public AbsKalmanFitter {

 private:

  DAF(const DAF&);
  DAF& operator=(genfit::DAF const&);

 public:

  /**
   * @brief Create DAF. Per default, use KalmanFitterRefTrack as fitter.
   *
   * @param useRefKalman If false, use KalmanFitter as fitter.
   */
  DAF(bool useRefKalman = true, double deltaPval = 1e-3, double deltaWeight = 1e-3);
  /**
   * @brief Create DAF. Use the provided AbsKalmanFitter as fitter.
   */
  DAF(AbsKalmanFitter* kalman, double deltaPval = 1e-3, double deltaWeight = 1e-3);
  ~DAF() {};

  //! Process a track using the DAF.
  void processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits = false);

  /** @brief Set the probability cut for the weight calculation for the hits.
   *
   * By default the cut values for measurements of dimensionality from 1 to 5 are calculated.
   * If you what to have cut values for an arbitrary measurement dimensionality use
   * addProbCut(double prob_cut, int maxDim);
   */
  void setProbCut(const double prob_cut);

  //! Set the probability cut for the weight calculation for the hits for a specific measurement dimensionality.
  void addProbCut(const double prob_cut, const int measDim);

  const std::vector<double>& getBetas() const {return betas_;}

  /** @brief Configure the annealing scheme.
   *
   * Set a start and end temperature and the number of steps. A logarithmic sequence of temperatures will be calculated.
   * Also sets #minIterations_ and #maxIterations_.
   */
  void setAnnealingScheme(double bStart, double bFinal, unsigned int nSteps);

  void setMaxIterations(unsigned int n) {maxIterations_ = n; betas_.resize(maxIterations_,betas_.back());}

  //! If all weights change less than delta between two iterations, the fit is regarded as converged.
  void setConvergenceDeltaWeight(double delta) {deltaWeight_ = delta;}

  AbsKalmanFitter* getKalman() const {return kalman_.get();}

  virtual void setMaxFailedHits(int val) {getKalman()->setMaxFailedHits(val);}

  virtual void setDebugLvl(unsigned int lvl = 1) {AbsFitter::setDebugLvl(lvl); if (lvl > 1) getKalman()->setDebugLvl(lvl-1);}

 private:

  /** @brief Calculate and set the weights for the next fitting pass.
   * Return if convergence is met.
   * The convergence criterium is the largest absolute change of all weights.
    */
  bool calcWeights(Track* trk, const AbsTrackRep* rep, double beta);


  double deltaWeight_; // convergence criterium
  std::vector<double> betas_;   // Temperatures, NOT inverse temperatures.
  double chi2Cuts_[7];  // '7' assumes tracks are helices with one
			// parameter, i.e. we're living in 3D space,
			// where time may be used in the fit.  Zeroth
			// entry is not used.

  std::unique_ptr<AbsKalmanFitter> kalman_;

 public:

  ClassDef(DAF,2)

};

}  /* End of namespace genfit */
/** @} */

#endif //genfit_DAF_h
