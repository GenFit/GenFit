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

#ifndef genfit_ThinScatterer_h
#define genfit_ThinScatterer_h

#include "SharedPlanePtr.h"
#include "Material.h"

#include <TObject.h>


namespace genfit {

/**
 * @brief Thin or thick scatterer
 */
class ThinScatterer : public TObject {

 public:

  ThinScatterer() :
    TObject(), sharedPlane_(), material_() {;}
  ThinScatterer(const SharedPlanePtr& sharedPlane, const Material& material) :
    TObject(), sharedPlane_(sharedPlane), material_(material) {;}

  SharedPlanePtr getPlane() const {return sharedPlane_;}
  const Material& getMaterial() const {return material_;}

  void Print(const Option_t* = "") const;

 private:

  SharedPlanePtr sharedPlane_; //! Material boundary.  '!' shuts up ROOT.
  Material material_; // Material properties


 public:
  ClassDef(ThinScatterer, 2)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_ThinScatterer_h
