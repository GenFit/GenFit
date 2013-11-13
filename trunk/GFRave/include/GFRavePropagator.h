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

/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup GFRave
 * @{
 */

#ifndef GFRAVEPROPAGATOR_H
#define GFRAVEPROPAGATOR_H

#include "GFRaveVertexFactory.h"
#include "AbsTrackRep.h"

#include <rave/Propagator.h>
#include <rave/Track.h>
#include <rave/Plane.h>
#include <rave/Cylinder.h>

#include <map>


namespace genfit {

/**
 * @brief GFRavePropagator class
 *
 * Inherits from rave::Propagator. A map has to be provided,
 * containing pointers the genfit::Tracks, as well as pointers to clones of fitted states.
 * The GFRavePropagator uses the information of the rave::track to set
 * the state and covariance of the corresponding MeasuredStateOnPlane,
 * extrapolates and then returns a new rave::track with the
 * extrapolated state and covariance.
 */
class GFRavePropagator : public rave::Propagator
{
  public:
    GFRavePropagator();
    virtual GFRavePropagator* copy() const;
    virtual rave::Track closestTo ( const rave::Track &,
                                    const rave::Point3D &, bool transverse ) const;
    virtual std::pair < rave::Track, double > to ( const rave::Track & orig,
                                                   const ravesurf::Plane & ) const;
    virtual std::pair < rave::Track, double > to ( const rave::Track & orig,
                                                   const ravesurf::Cylinder & ) const;

    virtual ~GFRavePropagator();

    void setIdGFTrackStateMap(std::map < int, genfit::trackAndState > * map);

  private:

    // data members
    std::map < int, genfit::trackAndState > * IdGFTrackStateMap_; // pointers to genfit::tracks and measuredStateOnPlanes via rave track ID
};

} /* End of namespace genfit */
/** @} */

#endif // GFRAVEPROPAGATOR_H
