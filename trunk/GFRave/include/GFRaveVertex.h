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

#ifndef GFRAVEVERTEX_H
#define GFRAVEVERTEX_H

#include "AbsTrackRep.h"
#include "Track.h"
#include "GFRaveTrackParameters.h"


namespace genfit {

/**
 * @brief GFRaveVertex class
 *
 * A Vertex contains information about its position and covariance.
 * The tracks the vertex is consisting of are stored in smoothedTracks_.
 * These GFRaveTrackParameters contain the weight of the corresponding track
 * in the vertex, smoothed track parameters and a pointer to the original
 * unaltered genfit::Track.
 */
class GFRaveVertex : public TObject {

  public:
    // constructors, destructors
    GFRaveVertex();
    GFRaveVertex(const TVector3 & pos, const TMatrixDSym & cov,
                 const std::vector < genfit::GFRaveTrackParameters* > & smoothedTracks,
                 double ndf, double chi2, int id = -1);

    GFRaveVertex(const GFRaveVertex &);

    GFRaveVertex& operator=(GFRaveVertex);
    void swap(GFRaveVertex&);

    ~GFRaveVertex();


    // Accessors
    //! get Position
    TVector3 getPos() const {return pos_;}

    //!get 3x3 covariance (error) of position.
    TMatrixDSym getCov() const {return cov_;}

    double getNdf() const {return ndf_;}
    double getChi2() const {return chi2_;}

    //! Number of tracks the vertex is made of
    unsigned int getNTracks() const {return smoothedTracks_.size();}
    GFRaveTrackParameters* getParameters(unsigned int i) const {return smoothedTracks_[i];}

    int getId() const {return id_;}

    void Print(const Option_t* = "") const;


  private:

    TVector3 pos_; // position of the vertex
    TMatrixDSym cov_; // error of the vertex position
    double ndf_;
    double chi2_;
    int id_; // id of the rave::vertex the GFVertex is created from

    std::vector < genfit::GFRaveTrackParameters* > smoothedTracks_; //-> track parameters of smoothed (with the vertex information) tracks, weights and original tracks; Vertex has ownership!

  public:
    ClassDef(GFRaveVertex, 1)

};

} /* End of namespace genfit */
/** @} */

#endif // GFRAVEVERTEX_H
