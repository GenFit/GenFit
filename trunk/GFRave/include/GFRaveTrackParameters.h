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

#ifndef GFRAVETRACKPARAMETERS_H
#define GFRAVETRACKPARAMETERS_H

#include "Track.h"
#include "AbsTrackRep.h"

#include <TObject.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TRef.h>

#include <iostream>


namespace genfit {

/**
 * @brief GFRaveTrackParameters class
 * Contains a pointer to the original genfit::Track, the weight of the track in the vertex,
 * and smoothed (with the vertex information) state and covariance of the track.
 */
class GFRaveTrackParameters : public TObject
{
  public:
    // constructors, destructors
    GFRaveTrackParameters();
    GFRaveTrackParameters(const Track* track, MeasuredStateOnPlane* originalState, double weight, const TVectorD & state6, const TMatrixDSym & cov6x6, bool isSmoothed);
    GFRaveTrackParameters(const Track* track, MeasuredStateOnPlane* originalState, double weight);

    // Accessors
    double getWeight() const {return weight_;}

    bool hasTrack() const {return originalTrack_.GetObject() != NULL;}
    const Track* getTrack() const {return  static_cast<Track*>(originalTrack_.GetObject());}

    UInt_t GetUniqueID() const {return originalTrack_.GetUniqueID();}

    bool hasSmoothedData() const {return hasSmoothedData_;}
    TVectorD getState() const {return state_;}
    TVector3 getPos() const;
    TVector3 getMom() const;
    const TMatrixDSym & getCov() const {return cov_;}

    double getCharge() const;
    double getPdg() const;

    void Print(const Option_t* = "") const;

  private:

    const TRef originalTrack_; // NO ownership. We use TRef, since the Tracks could be stored in another file or tree.

    double weight_; // weight of the track in the vertex
    TVectorD state_; // x, y, z, px, py, pz
    TMatrixDSym cov_; // 6x6 covariance matrix
    bool hasSmoothedData_; // true if state_ is forced to go through the vertex

  private:
    ClassDef(GFRaveTrackParameters, 1)
};

} /* End of namespace genfit */
/** @} */

#endif // GFRAVETRACKPARAMETERS_H
