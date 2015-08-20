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
    
    virtual GblTrackSegmentController* clone() const {return new GblTrackSegmentController(*this);}
    
    virtual void controlTrackSegment(TVector3, TVector3, GblFitter *) {;}
    
    virtual void Print(const Option_t* = "") const {;}
    
  protected:
    
  public:
    
    ClassDef(GblTrackSegmentController, 1)
    
  };
  
} /* End of namespace genfit */
/** @} */

#endif // genfit_GblTrackSegmentController_h
