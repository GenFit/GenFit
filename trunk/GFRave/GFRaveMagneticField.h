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

/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup GFRave
 * @{
 */

#ifndef GFRAVEMAGNETICFIELD_H
#define GFRAVEMAGNETICFIELD_H

// overwrite visual c stuff
#define RaveDllExport

#include <rave/MagneticField.h>

namespace gfrave
{

/** 
 * @brief GFRaveMagneticField class
 */

class GFRaveMagneticField : public rave::MagneticField
{
  public:
    GFRaveMagneticField(){};
    GFRaveMagneticField * copy() const;

    rave::Vector3D inTesla ( const rave::Point3D & ) const;

    //virtual ~GFRaveMagneticField();

  private:
};

}

#endif
