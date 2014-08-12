/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

#ifndef genfit_KalmanFitStatus_h
#define genfit_KalmanFitStatus_h

#include "FitStatus.h"

#include <Math/ProbFunc.h>


namespace genfit {

/**
 * @brief FitStatus for use with AbsKalmanFitter implementations
 */
class KalmanFitStatus : public FitStatus {

 public:

  KalmanFitStatus() :
    FitStatus(), numIterations_(0), fittedWithDaf_(false), fittedWithReferenceTrack_(false),
    trackLen_(0), fChi2_(-1e99), fNdf_(-1e99) {;}

  virtual ~KalmanFitStatus() {};

  virtual FitStatus* clone() const {return new KalmanFitStatus(*this);}

  unsigned int getNumIterations() const {return numIterations_;}
  bool isFittedWithDaf() const {return fittedWithDaf_;}
  bool isFittedWithReferenceTrack() const {return fittedWithReferenceTrack_;}
  double getTrackLen() const {return trackLen_;}
  double getForwardChi2() const {return fChi2_;}
  double getBackwardChi2() const {return FitStatus::getChi2();}
  double getForwardNdf() const {return fNdf_;}
  double getBackwardNdf() const {return FitStatus::getNdf();}
  // virtual double getPVal() : not overridden, as it does the right thing.
  double getForwardPVal() const {return std::max(0.,ROOT::Math::chisquared_cdf_c(fChi2_, fNdf_));}
  double getBackwardPVal() const {return FitStatus::getPVal(); }

  void setNumIterations(unsigned int numIterations) {numIterations_ = numIterations;}
  void setIsFittedWithDaf(bool fittedWithDaf = true) {fittedWithDaf_ = fittedWithDaf;}
  void setIsFittedWithReferenceTrack(bool fittedWithReferenceTrack = true) {fittedWithReferenceTrack_ = fittedWithReferenceTrack;}
  void setTrackLen(double trackLen) {trackLen_ = trackLen;}
  void setForwardChi2(double fChi2) {fChi2_ = fChi2;}
  void setBackwardChi2(double bChi2) {FitStatus::setChi2(bChi2);}
  void setForwardNdf(double fNdf) {fNdf_ = fNdf;}
  void setBackwardNdf(double bNdf) {FitStatus::setNdf(bNdf);}

  virtual void Print(const Option_t* = "") const;

 protected:

  unsigned int numIterations_; // number of iterations that have been performed
  bool fittedWithDaf_;
  bool fittedWithReferenceTrack_;

  double trackLen_;

  double fChi2_; // chi^2 of the forward fit
  double fNdf_; // degrees of freedom of the forward fit
  double fPval_; // p-value of the forward fit, set whenever either of chi2 or ndf changes

 public:

  ClassDef(KalmanFitStatus, 1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFitStatus_h
