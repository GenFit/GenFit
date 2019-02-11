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

#ifndef genfit_GblTrackSegmentController_h
#define genfit_GblTrackSegmentController_h

#include "Track.h"
#include "GblFitterInfo.h"

#include <Math/ProbFunc.h>
#include "GblTrajectory.h"
#include "GblFitter.h"
#include <TVector3.h>

namespace genfit {
  
  class GblFitter;
  /**
   * @brief TrackSegmentController for use with GblFitter
   */
  class GblTrackSegmentController {
    
  public:
    
    GblTrackSegmentController() {;}
    
    virtual ~GblTrackSegmentController() {};
    
    /**
    * @brief Function called for each segment of trajectory. User can decide on MS options.
    *        This function must be implemented by the actual class deriving from this abstract class
    * @param entry Position of segment starting point
    * @param exit Position of segment ending point
    * @param scatTheta Total MS variance accumulated in this segment
    * @param fitter Pointer to the fitter - so you can set the MS options
    */
    virtual void controlTrackSegment(TVector3 entry, TVector3 exit, double scatTheta, GblFitter * fitter) = 0;
    
    virtual void Print(const Option_t* = "") const {;}
    
  protected:
    
  public:
    
    ClassDef(GblTrackSegmentController, 1)
    
  };
  
} /* End of namespace genfit */
/** @} */

#endif // genfit_GblTrackSegmentController_h
