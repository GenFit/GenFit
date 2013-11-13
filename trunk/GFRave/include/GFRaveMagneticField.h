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

/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 */

/** @addtogroup GFRave
 * @{
 */

#ifndef GFRAVEMAGNETICFIELD_H
#define GFRAVEMAGNETICFIELD_H

#include <rave/MagneticField.h>


namespace genfit {

/**
 * @brief GFRaveMagneticField class
 * Uses the FieldManager to provide a magnetic field to rave.
 */
class GFRaveMagneticField : public rave::MagneticField {
  public:
    GFRaveMagneticField(){};
    virtual GFRaveMagneticField * copy() const;
    virtual ~GFRaveMagneticField(){};

    virtual rave::Vector3D inTesla ( const rave::Point3D & ) const;

  private:

};

} /* End of namespace genfit */
/** @} */

#endif // GFRAVEMAGNETICFIELD_H
