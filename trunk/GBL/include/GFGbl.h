/* Copyright 2013
   Authors: Sergey Yashchenko and Tadeas Bilka

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
/** @addtogroup genfit
 * @{
 */

#ifndef GFGBL_H
#define GFGBL_H

#include "GblTrajectory.h"
#include "AbsFitter.h"

#include <map>
#include <iostream>

#include <TMatrixD.h>
#include <assert.h>
#include <sstream>

#include <TMath.h>
#include <TVector3.h>


namespace genfit {


/** @brief Generic GBL implementation
 *
 * The interface class to GBL track fit
 *
 */
class GFGbl : public AbsFitter {

 private:
  GFGbl(const GFGbl&);
  GFGbl& operator=(GFGbl const&);

 public:

  /**
   * Constructor
   */
  GFGbl();

  /**
   * Destructor
   */
  virtual ~GFGbl() {;}

  /**
   * Performs fit on a Track.
   * Hit resorting currently NOT supported.
   */
  void processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits = false);


 public:

  ClassDef(GFGbl, 1)

};

}  /* End of namespace genfit */
/** @} */

#endif // GFGBL_H
