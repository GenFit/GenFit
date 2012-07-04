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

#include "GFRaveMagneticField.h"
#include "GFRavePropagator.h"

#include "GFException.h"

#include "rave/Vertex.h"
#include "rave/Ellipsoid3D.h"



GFRaveVertexFactory::GFRaveVertexFactory(int verbosity, bool useVacuumPropagator)
{
  fIdGFTrackMap = new std::map<int, GFTrack*>;
  fIdGFTrackRepMap = new std::map<int, GFAbsTrackRep*>;

  if (useVacuumPropagator) {
    fPropagator = new rave::VacuumPropagator();
  }
  else {
    fPropagator = new GFRavePropagator();
    ((GFRavePropagator*)fPropagator)->setIdGFTrackRepMap(fIdGFTrackRepMap);
  }

  fMagneticField = new GFRaveMagneticField();

  if (verbosity > 0) ++verbosity; // verbosity has to be >1 for rave

  fFactory = new rave::VertexFactory(*fMagneticField, *fPropagator, "default", verbosity); // here copies of fMagneticField and fPropagator are made!
}


GFRaveVertexFactory::~GFRaveVertexFactory(){
  clearMaps();
  delete fIdGFTrackMap;
  delete fIdGFTrackRepMap;

  delete fMagneticField;
  delete fPropagator;
  delete fFactory;
}


void
GFRaveVertexFactory::findVertices ( std::vector <  GFRaveVertex* > * GFvertices, const std::vector < GFTrack* > & GFTracks, bool use_beamspot ){
  clearMaps();

  try{
    GFRave::RaveToGFVertices(GFvertices,
                             fFactory->create(GFRave::GFTracksToTracks(GFTracks, fIdGFTrackMap, fIdGFTrackRepMap, 0),
                                              use_beamspot),
                             fIdGFTrackMap, fIdGFTrackRepMap);
  }
  catch(GFException & e){
    std::cerr << e.what();
    return;
  }
}


void
GFRaveVertexFactory::findVertices ( std::vector <  GFRaveVertex* > * GFvertices, const std::vector < GFAbsTrackRep* > & GFTrackReps, bool use_beamspot ){
  clearMaps();

  try{
    GFRave::RaveToGFVertices(GFvertices,
                             fFactory->create(GFRave::GFTrackRepsToTracks(GFTrackReps, fIdGFTrackMap, fIdGFTrackRepMap, 0),
                                              use_beamspot),
                             fIdGFTrackMap, fIdGFTrackRepMap);
  }
  catch(GFException & e){
    std::cerr << e.what();
    return;
  }
}


void
GFRaveVertexFactory::setBeamspot(const TVector3 & pos, const TMatrixT<double> & cov){
  fFactory->setBeamSpot(rave::Ellipsoid3D(GFRave::TVector3ToPoint3D(pos),
                        GFRave::TMatrixTToCovariance3D(cov)));
}


void
GFRaveVertexFactory::setMethod(const std::string & method){
  fFactory->setDefaultMethod(method);
  std::cout << "GFRaveVertexFactory::setMethod ==> set method to " << fFactory->method() << std::endl;
}


void
GFRaveVertexFactory::clearMaps(){
  fIdGFTrackMap->clear();

  for (unsigned int i=0; i<fIdGFTrackRepMap->size(); ++i){
    // in here are copies or trackreps -> we have ownership
    delete (*fIdGFTrackRepMap)[i];
  }
  fIdGFTrackRepMap->clear();
}
