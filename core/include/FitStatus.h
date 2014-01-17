/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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

#ifndef genfit_FitStatus_h
#define genfit_FitStatus_h

#include <Rtypes.h>
#include <Math/ProbFuncMathCore.h>


namespace genfit {

/** @brief Class where important numbers and properties of a fit can be stored.
 *
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 */
class FitStatus {

 public:

  FitStatus() :
    isFitted_(false), isFitConvergedFully_(false), isFitConvergedPartially_(false), nFailedPoints_(0),
    trackHasChanged_(false), trackIsPruned_(false), charge_(0), chi2_(-1e99), ndf_(-1e99)
  {;}

  virtual ~FitStatus() {};

  virtual FitStatus* clone() const {return new FitStatus(*this);}

  //! Has the track been fitted?
  bool isFitted() const {return isFitted_;}
  //! Did the fit converge (in all Points or only partially)?
  /**
   * Per default, this function will only be true, if all TrackPoints
   * (with measurements) have been used in the fit, and the fit has converged.
   *
   * If one or more TrackPoints have been skipped
   * (e.g. plane could not be constructed or extrapolation to plane failed),
   * but the fit otherwise met the convergence criteria, isFitConverged(false)
   * will return true.
   */
  bool isFitConverged(bool inAllPoints = true) const {
    if (inAllPoints) return isFitConvergedFully_; return isFitConvergedPartially_;
  }
  bool isFitConvergedFully() const {return isFitConvergedFully_;}
  bool isFitConvergedPartially() const {return isFitConvergedPartially_;}
  int getNFailedPoints() const {return nFailedPoints_;}
  //! Has anything in the Track been changed since the fit?
  bool hasTrackChanged() const {return trackHasChanged_;}
  //! Has the track been pruned after the fit?
  bool isTrackPruned() const {return trackIsPruned_;}
  //! Get the fitted charge.
  double getCharge() const {return charge_;}
  //! Get chi^2 of the fit.
  double getChi2() const {return chi2_;}
  //! Get the degrees of freedom of the fit.
  double getNdf() const {return ndf_;}
  /**
   * @brief Get the p value of the fit.
   *
   * Virtual, because the fitter may use a different probability distribution.
   */
  virtual double getPVal() const {return ROOT::Math::chisquared_cdf_c(chi2_, ndf_);}

  void setIsFitted(bool fitted = true) {isFitted_ = fitted;}
  void setIsFitConvergedFully(bool fitConverged = true) {isFitConvergedFully_ = fitConverged;}
  void setIsFitConvergedPartially(bool fitConverged = true) {isFitConvergedPartially_ = fitConverged;}
  void setNFailedPoints(int nFailedPoints) {nFailedPoints_ = nFailedPoints;}
  void setHasTrackChanged(bool trackChanged = true) {trackHasChanged_ = trackChanged;}
  void setIsTrackPruned(bool pruned = true) {trackIsPruned_ = pruned;}
  void setCharge(double charge) {charge_ = charge;}

  void setChi2(const double& chi2) {chi2_ = chi2;}
  void setNdf(const double& ndf) {ndf_ = ndf;}

  virtual void Print(const Option_t* = "") const;

 protected:

  //! has the track been fitted?
  bool isFitted_;
  //! did the fit converge with all TrackPoints?
  bool isFitConvergedFully_;
  //! did the fit converge with a subset of all TrackPoints?
  bool isFitConvergedPartially_;
  //! Number of failed TrackPoints
  int nFailedPoints_;
  //! has anything in the Track been changed since the fit? -> fit isn't valid anymore
  bool trackHasChanged_;
  //! Information has been stripped off, no refitting possible!
  bool trackIsPruned_;
  //! fitted charge
  double charge_;

  //! These are provided for the sake of the fitter, and their interpretation may vary.
  //! For the Kalman-derived fitters in particular, this corresponds to the backwards fit.
  double chi2_;
  double ndf_;
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_FitStatus_h
