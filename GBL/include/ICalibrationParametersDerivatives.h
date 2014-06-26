/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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

#ifndef genfit_ICalibrationParametersDerivatives_h
#define genfit_ICalibrationParametersDerivatives_h

#include "AbsMeasurement.h"
#include "StateOnPlane.h"
#include "TMatrixD.h"
#include <TMatrixT.h>


namespace genfit {

/** @brief Abstract base class to establish an interface between physical representation
 * of the detector for alignment/calibration and (fitted) state on genfit::Track
 * 
 * An implementing class (RecoHit) should be able to calculate derivatives of 
 * track fit residuals (2D for planar, 1D for strip hit) residuals w.r.t.
 * global (and optionally additional local) parameters. Global parameters
 * need unique integer labels to identify them across all sub-detectors
 *
 *  @author TadeasBilka
 */
class ICalibrationParametersDerivatives {

 public:
   virtual ~ICalibrationParametersDerivatives(){}

   /**
    * @brief Vector of integer labels for calibration/alignment
    * parameters available (must match #columns of derivatives(...))
    * 
    * unique across all sub-detectors in calibration
    * 
    * @return std::vector< int > Vector of integer labels
    */
   virtual std::vector<int> labels() = 0;
   
   /**
    * @brief Derivatives of residuals (local measurement coordinates) w.r.t. alignment/calibration parameters
    * Matrix "G" of derivatives valid for given prediction of track state:
    * 
    * G(i, j) = d_residual_i/d_parameter_j
    * 
    * For 2D measurement (u,v):
    * 
    * G = ( du/da du/db du/dc ... )
    *     ( dv/da dv/db dv/dc ... )
    * 
    * for calibration parameters a, b, c.
    * 
    * For 1D measurement both forms are allowed:
    * 
    * G = (   0     0     0   ... )
    *     ( dv/da dv/db dv/dc ... )    for V-strip,
    * 
    * 
    * G = ( du/da du/db du/dc ... )
    *     (   0     0     0   ... )    for U-strip,
    * 
    * or :
    * 
    * G = ( d_sensitive/da d_sensitive/db d_sensitive/dc ... )   as matrix with one row.
    * 
    * A possible algorithm using these derivatives
    * should be able to resolve this based on the measurement HMatrix.
    * Measurements with more dimesions (slopes, curvature) should provide
    * full 4-5Dx(n params) matrix (state as (q/p, u', v', u, v) or (u', v', u, v))
    * 
    * 
    * @param sop Predicted state of the track as linearization point around 
    * which derivatives of alignment/calibration parameters shall be computed
    * @return TMatrixD Matrix with #rows = dimension of residual, #columns = number of parameters.
    * #columns must match labels().size().
    */
   virtual TMatrixD derivatives(const genfit::StateOnPlane* sop) = 0;
   
   /**
    * @brief Derivatives for additional local parameters to be fitted
    * in global calibration algorithms together with with global parameters
    * 
    * Local parameters are not neccesarily identified by label because their number
    * is proportional to number of measurements included in calibration
    * (possibly very huge number!)
    * 
    * @return TMatrixD Matrix in form d_residual_i/d_parameter_j
    */
   virtual TMatrixD localDerivatives(const genfit::StateOnPlane*) {return TMatrixD();}
   
   /**
    * @brief Vector of integer labels for local calibration
    * parameters available (must match #columns of localDerivatives(...))
    * 
    * This will be usually ignored (e.g. does not have to match localDerivatives),
    * but it is a good practice to return vector of zeros of correct size
    * @return std::vector< int > Vector of integer labels
    */
   virtual std::vector<int> localLabels() {return std::vector<int>();}
   
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_ICalibrationParametersDerivatives_h
