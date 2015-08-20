/* Copyright 2008-2010, Technische Universitaet Muenchen,
 *   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch
 * 
 *   This file is part of GENFIT.
 * 
 *   GENFIT is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   GENFIT is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 * 
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_GblFitStatus_h
#define genfit_GblFitStatus_h

#include "FitStatus.h"
#include "Track.h"
#include "GblFitterInfo.h"

#include <Math/ProbFunc.h>
#include "GblTrajectory.h"


namespace genfit {
  
  /**
   * @brief FitStatus for use with GblFitter
   */
  class GblFitStatus : public FitStatus {
    
  public:
    
    GblFitStatus() :
    FitStatus(), numIterations_(0), fittedWithReferenceTrack_(false),
    trackLen_(0), curvatureFlag_(true), maxLocalFitParams_(0) {;}
    
    virtual ~GblFitStatus() {};
    
    virtual FitStatus* clone() const {return new GblFitStatus(*this);}
    
    void setCurvature(bool useCurvature) {curvatureFlag_ = useCurvature;}
    bool hasCurvature() { return curvatureFlag_; }
    void setMaxLocalFitParams(unsigned maxLocalFitParams) {maxLocalFitParams_ = maxLocalFitParams;}
    bool getMaxLocalFitParams() { return maxLocalFitParams_; }
    
    unsigned int getNumIterations() const {return numIterations_;}
    bool isFittedWithReferenceTrack() const {return fittedWithReferenceTrack_;}
    double getTrackLen() const {return trackLen_;}
    // virtual double getPVal() : not overridden, as it does the right thing.
    
    void setNumIterations(unsigned int numIterations) {numIterations_ = numIterations;}
    void setIsFittedWithReferenceTrack(bool fittedWithReferenceTrack = true) {fittedWithReferenceTrack_ = fittedWithReferenceTrack;}
    void setTrackLen(double trackLen) {trackLen_ = trackLen;}
    
    virtual void Print(const Option_t* = "") const {;}
    
  protected:
    
    unsigned int numIterations_; // number of iterations that have been performed
    bool fittedWithReferenceTrack_;
    
    double trackLen_;
    bool curvatureFlag_;
    int maxLocalFitParams_;
    
    
  public:
    
    ClassDef(GblFitStatus, 1)
    
  };
  
} /* End of namespace genfit */
/** @} */

#endif // genfit_GblFitStatus_h
