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


#ifndef GFPATHMAT_H
#define GFPATHMAT_H

#include "math.h"
#include "TVector3.h"
#include "TGeoMaterial.h"


class GFPathMat {

 public:
  // Constructors/Destructors
  GFPathMat();
  GFPathMat(double path, TGeoMaterial* mat){fPath = path; fMat = mat;}
  ~GFPathMat();

  // Accessors
  double getPath() const {return fPath;}
  double getAbsPath() const {return fabs(fPath);}
  TGeoMaterial* getMat() const {return fMat;}

  // Modifiers
  void setPath(double path){fPath = path;}
  void setMat(TGeoMaterial* const mat){fMat = mat;}

  // Functions
  void addToPath(double dpath){fPath += dpath;}
  void Print();

 private:
  double fPath; // pathlength to next position (signed)
  TGeoMaterial* fMat; // pointer to material
};

#endif

/** @} */
