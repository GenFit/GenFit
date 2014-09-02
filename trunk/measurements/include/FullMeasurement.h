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

#ifndef genfit_FullMeasurement_h
#define genfit_FullMeasurement_h

#include "AbsMeasurement.h"
#include "AbsHMatrix.h"
#include "MeasurementOnPlane.h"


namespace genfit {

class AbsTrackRep;

/** @brief Measurement class implementing a measurement of all track parameters.
 *
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This class can e.g. be used, if the fitted track parameters measured in one subdetector should be
 * put into one "measurement".
 */
class FullMeasurement : public AbsMeasurement {

 public:
  FullMeasurement(int nDim = 5);
  FullMeasurement(const MeasuredStateOnPlane&, int detId = -1, int hitId = -1, TrackPoint* trackPoint = NULL);

  virtual ~FullMeasurement() {;}

  virtual AbsMeasurement* clone() const {return new FullMeasurement(*this);}

  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const;

  virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const;

  virtual const AbsHMatrix* constructHMatrix(const AbsTrackRep*) const;

 protected:
  SharedPlanePtr plane_;   //! This is persistent, but '!' makes ROOT shut up.

 public:

  ClassDef(FullMeasurement,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_FullMeasurement_h
