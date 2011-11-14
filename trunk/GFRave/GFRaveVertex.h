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

#ifndef GFRAVEVERTEX_H
#define GFRAVEVERTEX_H

#include "GFAbsTrackRep.h"
#include "GFTrack.h"
#include "GFRaveTrackParameters.h"

#include <iostream>


/**
 * @brief GFRaveVertex class
 */

class GFRaveVertex : public TObject
{
  public:
    // constructors, destructors
    GFRaveVertex();
    GFRaveVertex(TVector3 pos, TMatrixT<double> cov,
                 std::vector < std::pair < double, GFTrack* > > originalTracks,
                 std::vector < std::pair < double, GFRaveTrackParameters > > smoothedTracks,
                 double ndf, double chi2, int id = -1);

    ~GFRaveVertex(){};

    // Modifiers


    // Accessors
    std::vector < std::pair < double, GFRaveVertex > >  getWeightedComponents() const {return fComponents;}
    std::pair < double, GFRaveVertex >  getWeightedComponents(unsigned int i) const {return fComponents[i];}

    unsigned int getNTracks() const {return fOriginalTracks.size();}
    unsigned int getNSmoothedTracks() const {return fSmoothedTracks.size();}
    unsigned int getNComponents() const {return fComponents.size();}

    TVector3 getPos() const {return fPos;}
    TMatrixT<double> getCov() const {return fCov;}

    std::vector < std::pair < double, GFTrack* > > getWeightedTracks() const {return fOriginalTracks;}
    std::pair < double, GFTrack* > getWeightedTracks(unsigned int i) const {return fOriginalTracks[i];}

    std::vector < double > getTrackWeights() const;
    double getTrackWeights(unsigned int i) const {return fOriginalTracks[i].first;}

    std::vector < GFTrack* > getTracks() const;
    GFTrack* getTracks(unsigned int i) const {return fOriginalTracks[i].second;}

    std::vector < std::pair < double, GFRaveTrackParameters > > getWeightedSmoothedParameters() const {return fSmoothedTracks;}
    std::pair < double, GFRaveTrackParameters > getWeightedSmoothedParameters(unsigned int i) const {return fSmoothedTracks[i];}

    std::vector < double > getSmoothedParametersWeights() const;
    double getSmoothedParametersWeights(unsigned int i) const {return fSmoothedTracks[i].first;}

    std::vector < GFRaveTrackParameters > getSmoothedParameters() const;
    GFRaveTrackParameters getSmoothedParameters(unsigned int i) const {return fSmoothedTracks[i].second;}


    double getNdf() const {return fNdf;}
    double getChi2() const {return fChi2;}

    int getId() const {return fId;}


  private:

    std::vector < std::pair < double, GFRaveVertex > > fComponents; // The vertex components - only used in the Gaussian algorithm.

    TVector3 fPos; // position of the vertex
    TMatrixT<double> fCov; // error of the vertex position
    std::vector < std::pair < double, GFTrack* > > fOriginalTracks; // original unaltered GFTracks with weights
    std::vector < std::pair < double, GFRaveTrackParameters > > fSmoothedTracks; // track parameters of smoothed (with the vertex information) tracks
    double fNdf;
    double fChi2;
    int fId; // id of the rave::vertex the GFVertex is created from

  private:
    ClassDef(GFRaveVertex,1)
};

#endif

/** @} */


