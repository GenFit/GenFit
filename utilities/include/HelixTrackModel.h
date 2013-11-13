/* Copyright 2013, Technische Universitaet Muenchen,
   Authors: Johannes Rauch

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


/** @addtogroup utilities
 * @{
 */

#ifndef genfit_HelixTrackModel_h
#define genfit_HelixTrackModel_h

#include <TObject.h>
#include <TVector3.h>


namespace genfit {

/**
 * @brief Helix track model for testing purposes
 */
class HelixTrackModel : public TObject {

 public:

  // Constructors/Destructors ---------
  HelixTrackModel(const TVector3& pos, const TVector3& mom, double charge);

  TVector3 getPos(double tracklength) const;
  void getPosMom(double tracklength, TVector3& pos, TVector3& mom) const;
  void getPosDir(double tracklength, TVector3& pos, TVector3& dir) const {
    getPosMom(tracklength, pos, dir);
    dir.SetMag(1);
  }


 private:

  double sgn_;
  double mom_;
  double R_; // radius
  TVector3 center_;
  double alpha0_;
  double theta_;


 public:
  ClassDef(HelixTrackModel,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HelixTrackModel_h
