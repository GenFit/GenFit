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

#include "TObjArray.h"

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
                 std::vector < GFRaveTrackParameters* > smoothedTracks,
                 double ndf, double chi2, int id = -1);

    GFRaveVertex(const GFRaveVertex &);

    GFRaveVertex& operator=(const GFRaveVertex & vertex);

    ~GFRaveVertex();


    // Modifiers


    // Accessors
    TVector3 getPos() const {return fPos;}
    TMatrixT<double> getCov() const {return fCov;}

    double getNdf() const {return fNdf;}
    double getChi2() const {return fChi2;}

    unsigned int getNTracks() const {return fSmoothedTracks->GetEntriesFast();}
    GFRaveTrackParameters* getParameters(unsigned int i) const {return (GFRaveTrackParameters*)fSmoothedTracks->At(i);}

    int getId() const {return fId;}

    void Print(const Option_t* = "") const;


  private:

    TVector3 fPos; // position of the vertex
    TMatrixT<double> fCov; // error of the vertex position
    double fNdf;
    double fChi2;
    int fId; // id of the rave::vertex the GFVertex is created from

    TObjArray* fSmoothedTracks; //-> track parameters of smoothed (with the vertex information) tracks, weights and original tracks

  private:
    ClassDef(GFRaveVertex,1)
};

#endif

/** @} */


