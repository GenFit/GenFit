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

/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup GFRave
 * @{
 */

#ifndef GFRAVEPROPAGATOR_H
#define GFRAVEPROPAGATOR_H

// overwrite visual c shit
#define RaveDllExport

#include <rave/Propagator.h>
#include <rave/Track.h>
#include <rave/Plane.h>
#include <rave/Cylinder.h>

#include <GFAbsTrackRep.h>
#include <GFTrack.h>

namespace gfrave
{

  /**
   * @brief GFRavePropagator class
   */

  class GFRavePropagator : public rave::Propagator
  {
    public:
      GFRavePropagator();
      virtual gfrave::GFRavePropagator* copy() const;
      virtual rave::Track closestTo ( const rave::Track &,
                                      const rave::Point3D &, bool transverse ) const;
      virtual std::pair < rave::Track, double > to ( const rave::Track & orig,
                                                     const ravesurf::Plane & ) const;
      virtual std::pair < rave::Track, double > to ( const rave::Track & orig,
                                                     const ravesurf::Cylinder & ) const;

      virtual ~GFRavePropagator();

      void setIdGFTrackMap(std::map<unsigned int, GFTrack*>* map){IdGFTrackMap = map;}

    private:


      // check if everything is ok, otherwise throw GFException; return corresponding Cardinal Rep
      GFAbsTrackRep* checkTrack(const rave::Track & track) const;

      // data members
      std::map<unsigned int, GFTrack*>* IdGFTrackMap; // pointers to GFTracks via rave track ID
  };

}

#endif
