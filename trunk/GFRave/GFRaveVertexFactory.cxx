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


gfrave::GFRaveVertexFactory::GFRaveVertexFactory(rave::Ellipsoid3D * beamspot)
{
  fMagneticField = new gfrave::GFRaveMagneticField();
  fPropagator = new gfrave::GFRavePropagator();

  if (beamspot==NULL) fFactory = new rave::VertexFactory(*fMagneticField, *fPropagator);
  else fFactory = new rave::VertexFactory(*fMagneticField, *fPropagator, *beamspot);
}


std::vector < rave::Vertex >
gfrave::GFRaveVertexFactory::create ( const std::vector < GFTrack* > & GFTracks, bool use_beamspot ) const{

  std::map<unsigned int, GFTrack*>* IdGFTrackMap;
  std::vector<rave::Track> ravetracks = gfrave::GFTracksToTracks(GFTracks, IdGFTrackMap, 0);
  fPropagator->setIdGFTrackMap(IdGFTrackMap);

  return fFactory->create(ravetracks);
}
