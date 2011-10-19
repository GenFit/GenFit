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


#include "GFRaveVertexFactory.h"
#include "GFRaveConverters.h"
#include "GFRaveVertex.h"

#include "GFException.h"

#include "rave/Vertex.h"
#include "rave/Ellipsoid3D.h"


GFRaveVertexFactory::GFRaveVertexFactory()
{
  fMagneticField = new GFRaveMagneticField();
  fPropagator = new GFRavePropagator();

  fFactory = new rave::VertexFactory(*fMagneticField, *fPropagator);
}


std::vector < GFRaveVertex* > *
GFRaveVertexFactory::create ( const std::vector < GFTrack* > & GFTracks, bool use_beamspot ) const{

  std::map<int, GFTrack*>* IdGFTrackMap; // bookkeeping of original GFTracks for later assignment to GFVertices
  std::map<int, GFAbsTrackRep*>* IdGFTrackRepMap; // map of copies of the cardinal reps for the GFRavePropagator; ownership of trackrep clones is HERE!!!

  std::vector < rave::Vertex > ravevertices;

  try{
    std::vector<rave::Track> ravetracks = GFRave::GFTracksToTracks(GFTracks, IdGFTrackMap, IdGFTrackRepMap, 0);
    fPropagator->setIdGFTrackMap(IdGFTrackRepMap);

    ravevertices = fFactory->create(ravetracks, use_beamspot);
  }
  catch(GFException & e){
    std::cerr << e.what();
  }

  return GFRave::RaveToGFVertices(ravevertices, IdGFTrackMap);
}


void
GFRaveVertexFactory::setBeamspot(const TVector3 & pos, const TMatrixT<double> & cov){
  fFactory->setBeamSpot(rave::Ellipsoid3D(GFRave::TVector3ToPoint3D(pos),
                        GFRave::TMatrixTToCovariance3D(cov)));
}


void
GFRaveVertexFactory::setMethod(const std::string & method){
  fFactory->setDefaultMethod(method);
}
