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

#include <TVector3.h>
#include <math.h>

class GFPointPath {

 public:
  // Constructors/Destructors
  GFPointPath();
  GFPointPath(const TVector3& pos, double path);
  GFPointPath(double posX, double posY, double posZ,
              double path);

  // Accessors
  double X() const {return fX;}
  double Y() const {return fY;}
  double Z() const {return fZ;}
  double getPath() const {return fPath;}
  double getAbsPath() const {return fabs(fPath);}

  //! get straight line distance to another pointPath
  double getDist(const GFPointPath& to) const;

  // Modifiers
  void setPos(const TVector3& pos){ fX=pos.X();  fY=pos.Y();  fZ=pos.Z();}
  void setPath(double path){fPath = path;}

  // Functions
  void addToPath(double dpath){fPath += dpath;}
  void Print() const;

 private:
  double fX,  fY,  fZ;  // position
  double fPath; // pathlength to next position (signed)
};

#endif

/** @} */
