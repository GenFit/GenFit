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


#ifndef GFRAVEVERTEXFACTORY_H
#define GFRAVEVERTEXFACTORY_H

#include "GFRaveVertex.h"
#include "Track.h"
#include "MeasuredStateOnPlane.h"

#include <vector>


namespace rave {
  class VertexFactory;
  class MagneticField;
  class Propagator;
}

namespace genfit {

/**
 * @brief Simple struct containing a Track pointer and a MeasuredStateOnPlane. Used in GFRave.
 */
struct trackAndState {
  const Track* track_; // pointer to original track
  MeasuredStateOnPlane* state_; // pointer to copy of fitted state; can be altered (extrapolated) during vertexing
};


/**
 * @brief Vertex factory for producing GFRaveVertex objects from Track objects.
 *
 * The GFRaveVertexFactory is basically a wrapper around the rave::VertexFactory.
 * It takes care of initializing the rave::VertexFactory, building the necessary maps,
 * convert GENFIT to rave objects and vice versa.
 **/
class GFRaveVertexFactory {
 public:
  // constructors, destructors
  GFRaveVertexFactory(int verbosity = 0, bool useVacuumPropagator = false);
  ~GFRaveVertexFactory();

  // functions
  void findVertices ( std::vector <  genfit::GFRaveVertex* > *, const std::vector < genfit::Track* > &, bool use_beamspot=false );
  //! MeasuredStateOnPlanes will be used (instead of the tracks fitted states) to calculate the rave::Track parameters. takes ownership of MeasuredStateOnPlanes.
  void findVertices ( std::vector <  genfit::GFRaveVertex* > *, const std::vector < genfit::Track* > &, std::vector < genfit::MeasuredStateOnPlane* > &, bool use_beamspot=false );

  void setBeamspot(const TVector3 & pos, const TMatrixDSym & cov);

  /**
   * Set the reconstruction method. See http://projects.hepforge.org/rave/trac/wiki/RaveMethods
   * Smoothing has to be turned on! e.g. kalman-smoothing:1
   */
  void setMethod(const std::string & method);

 private:

  void clearMap();

  // data members
  std::map<int, genfit::trackAndState> IdGFTrackStateMap_; // map of copies of the cardinal MeasuredStateOnPlanes for the GFRavePropagator; ownership of MeasuredStateOnPlanes is HERE!!!
  rave::VertexFactory* factory_; // Ownership
  rave::MagneticField* magneticField_; // Ownership
  rave::Propagator* propagator_; // Ownership

};

} /* End of namespace genfit */
/** @} */

#endif // GFRAVEVERTEXFACTORY_H
