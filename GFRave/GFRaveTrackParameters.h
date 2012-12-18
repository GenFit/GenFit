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

#ifndef GFRAVETRACKPARAMETERS_H
#define GFRAVETRACKPARAMETERS_H

#include "TObject.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "GFTrack.h"
#include "GFAbsTrackRep.h"

#include <iostream>


/**
 * @brief GFRaveTrackParameters class
 * Contains a pointer to the original GFTrack, the weight of the track in the vertex,
 * and smoothed (with the vertex information) state and covariance of the track.
 */

class GFRaveTrackParameters : public TObject
{
  public:
    // constructors, destructors
    GFRaveTrackParameters();
    GFRaveTrackParameters(GFTrack* track, GFAbsTrackRep* rep, double weight, const TVectorD & state6, const TMatrixDSym & cov6x6, bool isSmoothed);
    GFRaveTrackParameters(GFTrack* track, GFAbsTrackRep* rep, double weight);

    // Accessors
    double getWeight() const {return fWeight;}

    bool hasTrack() const {return fOriginalTrack!=NULL;}
    GFTrack* getTrack() const {return  fOriginalTrack;}

    bool hasRep() const {return fOriginalRep!=NULL;}
    GFAbsTrackRep* getRep() const {return  fOriginalRep;}

    bool hasSmoothedData() const {return fHasSmoothedData;}
    TVectorD getState() const {return fState;}
    TVector3 getPos() const;
    TVector3 getMom() const;
    const TMatrixDSym & getCov() const {return fCov;}

    double getCharge() const;
    double getPdg() const;

    void Print(const Option_t* = "") const;

  private:

    GFTrack* fOriginalTrack; // NO ownership
    GFAbsTrackRep* fOriginalRep; // NO ownership
    double fWeight;
    TVectorD fState; // x, y, z, px, py, pz
    TMatrixDSym fCov; // 6x6 covariance matrix
    bool fHasSmoothedData;

  private:
    ClassDef(GFRaveTrackParameters, 4)
};


/** @} */




#endif

/** @} */


