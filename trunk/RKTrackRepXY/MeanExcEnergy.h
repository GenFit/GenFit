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
#ifndef MEANEXCENERGY_H
#define MEANEXCENERGY_H

#include <iostream>
#include <assert.h>

class TGeoMaterial;

class MeanExcEnergy{
 private:
  MeanExcEnergy(){};
  static const int NELEMENTS = 92;
  static const float vals[];

 public:
  static float get(int Z){
    assert(Z>0&&Z<=NELEMENTS);
    return vals[Z-1];
  }
  static float get(TGeoMaterial*);

};


#endif
