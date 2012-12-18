/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "energyLoss.h"
#include "math.h"
#include "assert.h"
#include <iostream>

float energyLoss(float beta,float charge,float density,float ZoverA,float MEE){
  //std::cout << "beta " << beta << " charge " << charge << " density " << density << " ZoverA " << ZoverA<< " " << MEE << std::endl;
  float retVal = 0.3071*ZoverA*density/(beta*beta)*charge*charge;
  //std::cout << retVal << std::endl;
  //std::cout << (log(beta*beta*0.510999*2./(1.E-6*MEE*(1-beta*beta)))-beta*beta) << std::endl;
  retVal *= (log(beta*beta*0.510999*2./(1.E-6*MEE*(1-beta*beta)))-beta*beta);
  assert(retVal>0);
  //in GeV/cm, hence 1.e-3
  return 1.E-3*retVal;
}
