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

#ifndef genfit_TrackCandHit_h
#define genfit_TrackCandHit_h

#include <TObject.h>


namespace genfit {

/**
 * @brief Hit object for use in TrackCand. Provides IDs and sorting parameters.
 */
class TrackCandHit : public TObject {
 public:

  // Constructors/Destructors ---------
  TrackCandHit(int detId   = -1,
               int hitId   = -1,
               int planeId = -1,
               double sortingParameter  =  0.);

  virtual ~TrackCandHit() {;}

  virtual TrackCandHit* clone() const {return new TrackCandHit(*this);}

  // Accessors
  int    getDetId() const {return detId_;}
  int    getHitId() const {return hitId_;}
  int    getPlaneId() const {return planeId_;}
  double getSortingParameter() const {return sortingParameter_;}

  // Modifiers
  void setSortingParameter(double sortingParameter) {sortingParameter_ = sortingParameter;}

  virtual void Print(Option_t* option = "") const;


  /** @brief Equality operator. Does not check sortingParameter.
   */
  friend bool operator== (const TrackCandHit& lhs, const TrackCandHit& rhs);
  friend bool operator!= (const TrackCandHit& lhs, const TrackCandHit& rhs) {
    return !(lhs == rhs);
  }

  /** @brief Compare sortingParameter, needed for sorting
   */
  friend bool operator< (const TrackCandHit& lhs, const TrackCandHit& rhs) {
    return (lhs.sortingParameter_ < rhs.sortingParameter_);
  }

 protected:

  //! protect from calling copy c'tor from outside the class. Use #clone() if you want a copy!
  TrackCandHit(const TrackCandHit& other) :
    TObject(other), detId_(other.detId_), hitId_(other.hitId_), planeId_(other.planeId_), sortingParameter_(other.sortingParameter_) {;}
  //! protect from calling assignment operator from outside the class. Use #clone() instead!
  TrackCandHit& operator=(const TrackCandHit&);


  // Data Members ------------
  int    detId_; // detId id is -1 per default
  int    hitId_; // hitId id is -1 per default
  int    planeId_; // planeId id is -1 per default
  double sortingParameter_; // sorting parameter


 public:

  ClassDef(TrackCandHit,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackCandHit_h
