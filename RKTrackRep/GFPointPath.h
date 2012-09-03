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

/*
 */

/** @addtogroup RKTrackRep
 * @{
 */


#ifndef GFPOINTPATH_H
#define GFPOINTPATH_H

#include "math.h"
#include "TVector3.h"
#include "TGeoMaterial.h"


class GFPointPath {

 public:
  // Constructors/Destructors
  GFPointPath();
  GFPointPath(const TVector3& pos, double path){fPos = pos; fPath = path;}

  // Accessors
  const TVector3& getPos() const {return fPos;}
  double X() const {return fPos.X();}
  double Y() const {return fPos.Y();}
  double Z() const {return fPos.Z();}
  double getPath() const {return fPath;}
  double getAbsPath() const {return fabs(fPath);}

  // Modifiers
  void setPos(TVector3& pos){fPos = pos;}
  void setPath(double path){fPath = path;}

  // Functions
  void addToPath(double dpath){fPath += dpath;}
  void Print() const;

 private:
  TVector3 fPos; // position
  double fPath; // pathlength to next position (signed)
};

#endif

/** @} */
