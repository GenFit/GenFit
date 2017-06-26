/* Copyright 2008-2010, Technische Universitaet Muenchen,
 *   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch
 * 
 *   This file is part of GENFIT.
 * 
 *   GENFIT is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   GENFIT is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 * 
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_GblFitterInfo_h
#define genfit_GblFitterInfo_h

#include "AbsFitterInfo.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"
#include "StateOnPlane.h"
#include "TrackPoint.h"
#include "ThinScatterer.h"
#include "StateOnPlane.h"
#include "Track.h"
#include "ICalibrationParametersDerivatives.h"
#include "GblPoint.h"
#include "GblTrajectory.h"
#include "GblFitter.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "AbsHMatrix.h"

#include <vector>
#include <memory>


namespace genfit {
  
  
  /**
   *  @brief Collects information needed and produced by a GblFitter/GBL and is specific to one AbsTrackRep of the Track.
   */
  class GblFitterInfo : public AbsFitterInfo {
    
  public:
    
    /**
     * @brief Constructor for ROOT I/O 
     */
    GblFitterInfo();
    
    /**
     * @brief Default (inherited) constructor
     * Should not be used or the reference should set immediately upon
     * construction (to set the plane).
     */
    GblFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);
    
    /**
     * @brief Default user constructor
     * 
     * @param trackPoint The point at track to attach fitter info.
     * @param rep The representation this fitter info belongs to
     * @param referenceState State from extrapolation to init predictions and plane
     */
    GblFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep, StateOnPlane& referenceState);
    
    /**
     * @brief (Initial) reset of fitter info
     * 
     * @param measurementDim Measurement dimesion (2)
     * @param repDim Representation dimesion (5)
     * @return void
     */
    void reset(unsigned int measurementDim = 2, unsigned int repDim = 5);
        
    /**
     * @brief Set the prediction and plane (from measurement if any)
     * You should use the user constructor instead.
     * 
     * Reference gets not updated internally in fitter info.
     * After updateFitResults(...), updates affect the +/- predictions directly.
     * 
     * Should be called only once (so rather use constructor).
     * Otherwise will rewrite fitted state by reference (and you have to update from fit again)
     * 
     * @param referenceState StateOnPlane from extrapolation to this point
     * @return void
     */
    void setReferenceState(StateOnPlane& referenceState);
    /**
     * @brief Set the Jacobian for further GblPoint construction
     * 
     * @param jacobian 5x5 TMatrixD with Jacobian for propagation of the state from previous point
     *                 to this point.
     * @return void
     */
    void setJacobian(TMatrixD jacobian);
    
    /**
     * @brief Get scattering covariance projected into (measurement) plane
     * 
     * @param variance Variance of slopes in track frame
     * @param trackDirection Direction of the track at the plane
     * @param measurementPlane The plane with measurement to which MS shall be projected
     * @return TMatrixDSym
     */
    TMatrixDSym getCovariance(double variance, TVector3 trackDirection, SharedPlanePtr measurementPlane) const;
    
    /**
     * @brief Collect all data and create a GblPoint
     * - Jacobian is taken from fitter info ... use setJaobian() to change it
     * - A scatterer is added if ThinScatterer at point
     * - A measurement (first MeasurementOnPlane in 1st RawMeasurement) is added.
     * - Global and local derivatives are added for RawMesurement implementing
     *   ICalibrationParametersDerivatives interface. Using most up to date prediction.
     * 
     * @return gbl::GblPoint
     */
    gbl::GblPoint constructGblPoint();
  
    /**
     * @brief Update fitter info from GBL fit results
     * 
     * Locates itself in track/trajectory and updates
     * predictions, errors, etc.
     * 
     * @param traj Fitted GblTrajectory constructed with this point
     * @return void
     */
    void updateFitResults(gbl::GblTrajectory& traj);
    
    /**
     * @brief Get the prediction at this point
     * Always biased in GBL (global fit)
     * There are 2 states, before and after kink (if there is
     * a scatterer at this point).
     * Per default the state after kink for forward propagation
     * is returned.
     * 
     * If not yet fitted, returns the reference state
     * 
     * @param afterKink If true, returns prediction for forward propagation.
     *                  If false, for backward
     * @return const genfit::MeasuredStateOnPlane&
     */
    const MeasuredStateOnPlane& getFittedState(bool afterKink = true) const;

    /**
     * @brief Get the residual
     * 
     * Temporarily constructs measurements
     * and calculates residual as meas - prediction.
     * 
     * Always biased. Always only one (1st) measurement!
     * 
     * @param  ...
     * @param  ...
     * @param onlyMeasurementErrors If true, covariance of measurement returned.
     *                              If false, diagonalized residual error incl. correlation from track fit is returned.
     * @return genfit::MeasurementOnPlane
     */
    MeasurementOnPlane getResidual(unsigned int = 0, bool = false, bool onlyMeasurementErrors = true) const;
    
    /**
     * @brief Get kink (residual) with diagonalized covariance (2D)
     * Covariance may be zero if not yet fitted or no scatterer!
     * 
     * @return genfit::MeasurementOnPlane
     */
    MeasurementOnPlane getKink() const;
    
    /**
     * @brief Get kink (residual) (2D)
     * = 0 - ( (+)pred - (-)pred )
     * 
     * @return TVectorD
     */
    TVectorD getKinks() const;
    
    /**
     * @brief Get the measurement on plane from stored
     * measurement data (from last construction/update)
     * 
     * @return genfit::MeasuredStateOnPlane
     */
    MeasurementOnPlane getMeasurement() const;
    
    /**
     * @brief SHOULD BE USED ONLY INTERNALY!
     * Update the plane from measurement constructed with state
     * or take plane from state if there is no measurement.
     * 1st raw Measurement with highest weight is constructed and stored as matrices.
     * 
     * @return void
     */
    void updateMeasurementAndPlane(const StateOnPlane & sop);

    /**
     * @brief Returns (copy of) the stored reference 5D state at current plane with internal rep
     * 
     * @return genfit::StateOnPlane
     */
    StateOnPlane getReferenceState() const {return StateOnPlane(refPrediction_, sharedPlane_, rep_);}
    
    /**
     * @brief Re-extrapolates between prevFitterInfo and this point using
     * forward state to update the Jacobian (if planes and/or states changed,
     * internal predictions are extrapolated to new planes)
     * 
     * @param prevFitterInfo Pointer to GblFitterInfo of previous point
     * @return void
     */
    void recalculateJacobian(GblFitterInfo* prevFitterInfo);
    
    virtual ~GblFitterInfo() {;}
    virtual GblFitterInfo* clone() const;
    bool hasMeasurements() const {return trackPoint_->hasRawMeasurements();}
    bool hasReferenceState() const {return (refPrediction_(0) != 0.);}
    bool hasForwardPrediction() const {return hasReferenceState();}
    bool hasBackwardPrediction() const {return hasReferenceState();}
    bool hasForwardUpdate() const {return hasForwardPrediction();}
    bool hasBackwardUpdate() const {return hasBackwardPrediction();}
    bool hasUpdate(int direction) const {if (direction < 0) return hasBackwardPrediction(); return hasForwardPrediction();}
    bool hasPredictionsAndUpdates() const {return (hasForwardPrediction() && hasBackwardPrediction() && hasForwardUpdate() && hasBackwardUpdate());}

    void deleteForwardInfo() {;}
    void deleteBackwardInfo() {;}
    void deletePredictions() {
      deleteBackwardInfo();
      deleteForwardInfo();
    }
    void deleteReferenceInfo() {;} // Empty because we really do not want to delete reference without a new one
    void deleteMeasurementInfo() {;} // We do not keep the measurements
    virtual void Print(const Option_t* = "") const;
    virtual bool checkConsistency(const genfit::PruneFlags* = nullptr) const;
       
  private:
    TMatrixD jacobian_;
    TVectorD measResiduals_;
    TVectorD measResidualErrors_;
    TVectorD kinkResiduals_;
    TVectorD kinkResidualErrors_;
    TVectorD measDownWeights_;
    TVectorD kinkDownWeights_;
    TVectorD bwdStateCorrection_;
    TVectorD fwdStateCorrection_;
    TMatrixDSym bwdCov_;
    TMatrixDSym fwdCov_;
    TVectorD fwdPrediction_;
    TVectorD bwdPrediction_;
    TVectorD refPrediction_;
    
    TVectorD measurement_;
    TMatrixDSym measCov_;
    TMatrixD hMatrix_;

    mutable std::unique_ptr<MeasuredStateOnPlane> fittedStateBwd_; //!  cache
    mutable std::unique_ptr<MeasuredStateOnPlane> fittedStateFwd_; //!  cache
    
  public:
    
    ClassDef(GblFitterInfo, 1)
    
  };
  
} /* End of namespace genfit */
/** @} */

#endif // genfit_GblFitterInfo_h
