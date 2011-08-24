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

#include<iostream>
#include<cassert>
#include"GFRectFinitePlane.h"

GFRectFinitePlane::GFRectFinitePlane(const double& umin,const double& umax,
				     const double& vmin,const double& vmax)
  : fUmin(umin),fUmax(umax),fVmin(vmin),fVmax(vmax)
{assert(umin<umax);assert(vmin<vmax);}

GFRectFinitePlane::GFRectFinitePlane()
  : fUmin(1.),fUmax(-1.),fVmin(1.),fVmax(-1.)//for this default ctor inActive always false
{}


GFRectFinitePlane::~GFRectFinitePlane(){
  
}

bool GFRectFinitePlane::inActive(const double& u,const double& v) const{
  if(u>=fUmin && u<=fUmax && v>=fVmin && v<=fVmax) return true;
  return false;
}

void GFRectFinitePlane::Print(const Option_t* option) const{
  std::cout << "Rectangular Finite Plane Umin=" << fUmin << ", Umax="
	    << fUmax << ", Vmin=" << fVmin << ", Vmax=" << fVmax << std::endl;
};
