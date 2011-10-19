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
#include "TMatrixT.h"
#include "TVector3.h"

#include <iostream>


/**
 * @brief GFRaveTrackParameters class
 * Contains state and covariance of the tracks, smoothed with the vertex information
 */

class GFRaveTrackParameters : public TObject
{
  public:
    // constructors, destructors
    GFRaveTrackParameters();
    GFRaveTrackParameters(TMatrixT<double> state6, TMatrixT<double> cov6x6, double charge, int pdg);

    // member functions

    // Accessors
    TMatrixT<double> getState() const {return fState;}
    TVector3 getPos() const;
    TVector3 getMom() const;

    const TMatrixT<double> & getCov() const {return fCov;}
    double getCharge() const {return fCharge;}
    double getPdg() const {return fPdg;}

  private:

    TMatrixT<double> fState; // x, y, z, px, py, pz
    TMatrixT<double> fCov; // 6x6 covariance matrix
    double fCharge;
    double fPdg;

  private:
    ClassDef(GFRaveTrackParameters, 1)
};


/** @} */




#endif

/** @} */


