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
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_WireTrackCandHit_h
#define genfit_WireTrackCandHit_h

#include "TrackCandHit.h"


namespace genfit {

/**
 * @brief Hit object for use in TrackCand. Provides additional left/right parameter.
 */
class WireTrackCandHit : public TrackCandHit {
 public:

  // Constructors/Destructors ---------
  WireTrackCandHit(int detId   = -1,
               int hitId   = -1,
               int planeId = -1,
               double sortingParameter  =  0.,
               char leftRight = 0);

  virtual ~WireTrackCandHit() {;}

  virtual WireTrackCandHit* clone() const {return new WireTrackCandHit(*this);}

  // Accessors
  int getLeftRightResolution() const {return leftRight_;}

  // Modifiers
  /**
   * select how to resolve the left/right ambiguity:
   * -1: negative (left) side on vector (track direction) x (wire direction)
   * 0: auto select (take side with smallest distance to track)
   * 1: positive (right) side on vector (track direction) x (wire direction)
   */
  void setLeftRightResolution(int leftRight){
    if (leftRight==0) leftRight_ = 0;
    else if (leftRight<0) leftRight_ = -1;
    else leftRight_ = 1;
  }

  virtual void Print(Option_t* option = "") const;


 protected:

  //! protect from calling copy c'tor from outside the class. Use #clone() if you want a copy!
  WireTrackCandHit(const WireTrackCandHit& other) :
    TrackCandHit(other), leftRight_(other.leftRight_) {;}
  //! protect from calling assignment operator from outside the class. Use #clone() instead!
  WireTrackCandHit& operator=(const WireTrackCandHit&);


  // Data Members ------------
  signed char leftRight_; // left/right ambiguity handling


 public:

  ClassDef(WireTrackCandHit, 2)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackCandHit_h
