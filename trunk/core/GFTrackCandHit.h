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
/** @addtogroup genfit
 * @{ */

#ifndef GFTRACKCANDHIT_H
#define GFTRACKCANDHIT_H

#include <TObject.h>

/** @brief Track candidate Hit. Collection of indices.
 *
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 */

class GFTrackCandHit : public TObject {
 public:

  virtual GFTrackCandHit* clone() const {return new GFTrackCandHit(*this);}

  // Constructors/Destructors ---------
  GFTrackCandHit(int detId   = -1,
                 int hitId   = -1,
                 int planeId = -1,
                 double rho  =  0.);

  ~GFTrackCandHit();

  /** @brief Equality operator. Does not check rho.
   */
  friend bool operator== (const GFTrackCandHit& lhs, const GFTrackCandHit& rhs);
  friend bool operator!= (const GFTrackCandHit& lhs, const GFTrackCandHit& rhs) {
    return !(lhs == rhs);
  }

  /** @brief Compare rho, needed for sorting
   */
  friend bool operator< (const GFTrackCandHit& lhs, const GFTrackCandHit& rhs) {
    return (lhs.fRho < rhs.fRho);
  }

  // Accessors
  int    getDetId() const {return fDetId;}
  int    getHitId() const {return fHitId;}
  int    getPlaneId() const {return fPlaneId;}
  double getRho() const {return fRho;}

  virtual void Print(Option_t* option = "") const;

  // Modifiers
  void setRho(double rho) {fRho = rho;}

 protected:
  // Private Data Members ------------
  int    fDetId; // detId id is -1 per default
  int    fHitId; // hitId id is -1 per default
  int    fPlaneId; // planeId id is -1 per default
  double fRho; // sorting parameter

 public:
  ClassDef(GFTrackCandHit, 1)
};

#endif

/** @} */
