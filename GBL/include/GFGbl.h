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
  
  std::string m_milleFileName;
  std::string m_gblInternalIterations;
  double m_pValueCut;
  int m_minNdf;
  

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
   * Creates the mille binary file for output of
   * data for Millepede II alignment, can be set by setMP2Options
   */
  void beginRun(); 
  
  /**
   * Required to write and close ROOT file
   * with debug output. Destructor cannot be used.
   * To be called from endRun function of a module
   */
  void endRun();

  void setGBLOptions(std::string internalIterations = "THC", unsigned int externalIterations = 1);


  void setMP2Options(double pValueCut = 0., unsigned int minNdf = 1, std::string mille_file_name = "millefile.dat");

  
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

