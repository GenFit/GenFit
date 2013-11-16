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

#include "HelixTrackModel.h"
#include <FieldManager.h>

#include <assert.h>
#include <math.h>

namespace genfit {

HelixTrackModel::HelixTrackModel(const TVector3& pos, const TVector3& mom, double charge) {

  mom_ = mom.Mag();

  TVector3 B = genfit::FieldManager::getInstance()->getFieldVal(pos);

  // B must point in Z direction
  assert(B.X() == 0);
  assert(B.Y() == 0);

  double Bz = B.Z();

  // calc helix parameters
  TVector3 dir2D(mom);
  dir2D.SetZ(0);
  dir2D.SetMag(1.);
  R_ = 100.*mom.Perp()/(0.0299792458*Bz) / fabs(charge);
  sgn_ = 1;
  if (charge<0) sgn_=-1.;
  center_ = pos + sgn_ * R_ * dir2D.Orthogonal();
  alpha0_ = (pos-center_).Phi();

  theta_ = mom.Theta();

  //std::cout<<"radius " << R_ << "  center ";
  //center_.Print();

}


TVector3 HelixTrackModel::getPos(double tracklength) const {

  TVector3 pos;

  double angle = alpha0_ - sgn_ * tracklength / R_ * sin(theta_);

  TVector3 radius(R_,0,0);
  radius.SetPhi(angle);
  pos = center_ + radius;
  pos.SetZ(center_.Z() - sgn_ * ((alpha0_-angle)*R_ * tan(theta_-M_PI/2.)) );

  return pos;
}

void HelixTrackModel::getPosMom(double tracklength, TVector3& pos, TVector3& mom) const {

  double angle = alpha0_ - sgn_ * tracklength / R_ * sin(theta_);

  TVector3 radius(R_,0,0);
  radius.SetPhi(angle);
  pos = center_ + radius;
  pos.SetZ(center_.Z() - sgn_ * ((alpha0_-angle)*R_ * tan(theta_-M_PI/2.)) );

  mom.SetXYZ(1,1,1);
  mom.SetTheta(theta_);
  mom.SetPhi(angle - sgn_*M_PI/2.);
  mom.SetMag(mom_);

  /*std::cout<<"tracklength " << tracklength << "\n";
  std::cout<<"angle " << angle << "\n";
  std::cout<<"radius vector "; radius.Print();
  std::cout<<"pos "; pos.Print();
  std::cout<<"mom "; mom.Print();*/

}


} /* End of namespace genfit */
