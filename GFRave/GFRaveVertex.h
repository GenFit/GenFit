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
                 std::vector < GFRaveTrackParameters > smoothedTracks,
                 double ndf, double chi2, int id = -1);

    ~GFRaveVertex(){};

    // Modifiers


    // Accessors
    //std::vector < std::pair < double, GFRaveVertex > >  getWeightedComponents() const {return fComponents;}
    //std::pair < double, GFRaveVertex >  getWeightedComponents(unsigned int i) const {return fComponents[i];}

    TVector3 getPos() const {return fPos;}
    TMatrixT<double> getCov() const {return fCov;}

    double getNdf() const {return fNdf;}
    double getChi2() const {return fChi2;}

    unsigned int getNTracks() const {return fSmoothedTracks.size();}
    std::vector < GFRaveTrackParameters > getParameters() const {return fSmoothedTracks;}
    GFRaveTrackParameters getParameters(unsigned int i) const {return fSmoothedTracks[i];}

    int getId() const {return fId;}

    void Print(const Option_t* = "") const;


  private:

    //std::vector < std::pair < double, GFRaveVertex > > fComponents; // The vertex components - only used in the Gaussian algorithm.

    TVector3 fPos; // position of the vertex
    TMatrixT<double> fCov; // error of the vertex position
    double fNdf;
    double fChi2;
    int fId; // id of the rave::vertex the GFVertex is created from

    std::vector < GFRaveTrackParameters > fSmoothedTracks; // track parameters of smoothed (with the vertex information) tracks, weights and original tracks

  private:
    ClassDef(GFRaveVertex,1)
};

#endif

/** @} */


