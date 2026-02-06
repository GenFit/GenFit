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

#ifndef genfit_GblFitterInfo2_h
#define genfit_GblFitterInfo2_h

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
#include "GblFitter2.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "AbsHMatrix.h"

#include <vector>
#include <memory>


namespace genfit {
  
  
  /**
   *  @brief Collects information needed and produced by a GblFitter/GBL and is specific to one AbsTrackRep of the Track.
   */
  class GblFitterInfo2 : public AbsFitterInfo {
    
  public:
    
    /**
     * @brief Constructor for ROOT I/O 
     */
    GblFitterInfo2();
    
    /**
     * @brief Default (inherited) constructor
     * Should not be used or the reference should set immediately upon
     * construction (to set the plane).
     */
    GblFitterInfo2(const TrackPoint* trackPoint, const AbsTrackRep* rep);
    
    /**
     * @brief Default user constructor
     * 
     * @param trackPoint The point at track to attach fitter info.
     * @param rep The representation this fitter info belongs to
     * @param referenceState State from extrapolation to init predictions and plane
     */
    GblFitterInfo2(const TrackPoint* trackPoint, const AbsTrackRep* rep, StateOnPlane& referenceState);
    
    /**
     * @brief (Initial) reset of fitter info
     * 
     * @param measurementDim Measurement dimesion (2)
     * @param repDim Representation dimesion (5)
     * @return void
     */
    void reset(unsigned int repDim = 5);
        
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
     * @brief Collect all data and create a GblPoint
     * - Jacobian is taken from fitter info ... use setJaobian() to change it
     * - A scatterer is added if ThinScatterer at point
     * - A measurement (first MeasurementOnPlane in 1st RawMeasurement) is added.
     * - Global and local derivatives are added for RawMesurement implementing
     *   ICalibrationParametersDerivatives interface. Using most up to date prediction.
     * 
     * @param index  measurement index
     * @param allowAmbiguities Allow ambiguities (to be resolved by GBL)
     * @return gbl::GblPoint
     */
    gbl::GblPoint constructGblPoint(unsigned int index, bool allowAmbiguities);
  
    /**
     * @brief Resolve ambuguities from GBL fit
     * 
     * @param traj Fitted GblTrajectory 
     * @param point GblPoint, optionally with ambiguties to be resolved
     * @param nResolved  number of resolved ambiguities
     * @param nSwapped   number of swaps
     * @return void
     */
    void resolveAmbiguities(gbl::GblTrajectory& traj, gbl::GblPoint& point, unsigned int &nResolved, unsigned int &nSwapped);

    /**
     * @brief Update down weights from GBL fit
     * 
     * @param traj Fitted GblTrajectory 
     * @return void
     */
    void updateDownweights(gbl::GblTrajectory& traj);

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
    const MeasuredStateOnPlane& getFittedState(bool afterKink = true) const override;

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
    MeasurementOnPlane getResidual(unsigned int = 0, bool = false, bool = false) const override;
    
    /**
     * @brief Get kinks and steps (residual) (4D)
     * = 0 - ( (+)pred - (-)pred )
     * 
     * @return TVectorD
     */
    TVectorD getKinksAndSteps() const;
    
    /**
     * @brief Get the measurement on plane from stored
     * measurement data (from last construction/update)
     *
     * @param Measurement (ambiguity)
     *
     * @return genfit::MeasuredStateOnPlane
     */
    MeasurementOnPlane getMeasurement(unsigned int = 0) const;
    
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
    void recalculateJacobian(GblFitterInfo2* prevFitterInfo);
    
    /**
     * @brief Set GBL label
     * 
     * @param unsigned int  label (>0)
     * @return void
     */
    void setLabel(unsigned int aLabel) {label_ = aLabel;}
    
    /**
     * @brief Get GBL label
     * 
     * @return unsigned int
     */    
    unsigned int getLabel() const {return label_;}

    /**
     * @brief Set GBL downweight
     * 
     * @param unsigned int  measurement index
     * @param double  downweight
     * @return void
     */
    void setDownWeight(unsigned int iMeas, double downWeight) {measDownWeight_.at(iMeas) = downWeight;}
    
    /**
     * @brief Get GBL downWeight
     * 
     * @param unsigned int  measurement index
     * @return double
     */    
    double getDownWeight(unsigned int iMeas) const {return measDownWeight_.at(iMeas);}

    /**
     * @brief Get GBL downWeight
     * 
     * @param unsigned int  measurement index
     * @return double
     */    
    const std::vector<double>& getDownWeights() const {return measDownWeight_;}
    
    /**
     * @brief Get most probable measurement (index)
     * 
     * @return unsigned int
     */    
    unsigned int getMopMeas() const {return measUsed_.at(0);}
    
    virtual ~GblFitterInfo2() {;}
    virtual GblFitterInfo2* clone() const override;
    bool hasMeasurements() const override {return trackPoint_->hasRawMeasurements();}
    bool hasReferenceState() const override {return (refPrediction_(0) != 0.);}
    bool hasForwardPrediction() const override {return hasReferenceState();}
    bool hasBackwardPrediction() const override {return hasReferenceState();}
    bool hasForwardUpdate() const override {return hasForwardPrediction();}
    bool hasBackwardUpdate() const override {return hasBackwardPrediction();}
    bool hasUpdate(int direction) const override {if (direction < 0) return hasBackwardPrediction(); return hasForwardPrediction();}
    bool hasPredictionsAndUpdates() const {return (hasForwardPrediction() && hasBackwardPrediction() && hasForwardUpdate() && hasBackwardUpdate());}

    void deleteForwardInfo() override {;}
    void deleteBackwardInfo() override {;}
    void deletePredictions() {
      deleteBackwardInfo();
      deleteForwardInfo();
    }
    void deleteReferenceInfo() override {;} // Empty because we really do not want to delete reference without a new one
    void deleteMeasurementInfo() override {;} // We do not keep the measurements
    virtual void Print(const Option_t* = "") const override;
    virtual bool checkConsistency(const genfit::PruneFlags* = nullptr) const override;
       
  private:
    enum measurementType{_uMeas, _vMeas, _uvMeas, _noMeas};
    measurementType measType_;
    unsigned int label_;
    TMatrixD jacobian_;
    TMatrixDSym noise_;
    TVectorD bwdStateCorrection_;
    TVectorD fwdStateCorrection_;
    TMatrixDSym bwdCov_;
    TMatrixDSym fwdCov_;
    TVectorD diffFwdBwdPrediction_;
    TVectorD fwdPrediction_;
    TVectorD bwdPrediction_;
    TVectorD refPrediction_;
    TMatrixD hMatrix_;

    // allow for left/right ambiguity in CDC
    std::vector<unsigned int> measUsed_; // list of used measurements (indices, sorted by prevoius weights)
    std::vector<double> measDownWeight_; // downweights from GBL (not normalized)
    std::vector<TVectorD> measurement_;
    std::vector<TMatrixDSym> measCov_;

    mutable std::unique_ptr<MeasuredStateOnPlane> fittedStateBwd_; //!  cache
    mutable std::unique_ptr<MeasuredStateOnPlane> fittedStateFwd_; //!  cache
    
  public:
    
    ClassDefOverride(GblFitterInfo2, 1)
    
  };
  
} /* End of namespace genfit */
/** @} */

#endif // genfit_GblFitterInfo2_h
