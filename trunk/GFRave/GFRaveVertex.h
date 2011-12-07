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
 * A Vertex contains information about its position and covariance.
 * The tracks the vertex is consisting of are stored in fSmoothedTracks.
 * These GFRaveTrackParameters contain the weight of the corresponding track
 * in the vertex, smoothed track parameters and a pointer to the original
 * unaltered GFTrack.
 */

class GFRaveVertex : public TObject
{
  public:
    // constructors, destructors
    GFRaveVertex();
    GFRaveVertex(const TVector3 & pos, const TMatrixT<double> & cov,
                 const std::vector < GFRaveTrackParameters* > & smoothedTracks,
                 double ndf, double chi2, int id = -1);

    GFRaveVertex(const GFRaveVertex &);

    GFRaveVertex& operator=(const GFRaveVertex &);

    ~GFRaveVertex();


    // Modifiers


    // Accessors
    /**
     * get Position
     */
    TVector3 getPos() const {return fPos;}

    /**
     * get 3x3 covariance (error) of position.
     */
    TMatrixT<double> getCov() const {return fCov;}

    double getNdf() const {return fNdf;}
    double getChi2() const {return fChi2;}

    /**
     * Number of tracks the vertex is made of
     */
    unsigned int getNTracks() const {return fSmoothedTracks.size();}
    GFRaveTrackParameters* getParameters(unsigned int i) const {return fSmoothedTracks[i];}

    int getId() const {return fId;}

    void Print(const Option_t* = "") const;


  private:

    TVector3 fPos; // position of the vertex
    TMatrixT<double> fCov; // error of the vertex position
    double fNdf;
    double fChi2;
    int fId; // id of the rave::vertex the GFVertex is created from

    std::vector < GFRaveTrackParameters* > fSmoothedTracks; // track parameters of smoothed (with the vertex information) tracks, weights and original tracks

  private:
    ClassDef(GFRaveVertex,1)
};

#endif

/** @} */


