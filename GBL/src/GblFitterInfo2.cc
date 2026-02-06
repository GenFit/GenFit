#include "GblFitterInfo2.h"
#include "MeasurementOnPlane.h"
#include "HMatrixUV.h"
#include "HMatrixV.h"
#include "HMatrixU.h"

//! Namespace for the general broken lines package
namespace genfit {

  GblFitterInfo2::GblFitterInfo2() : AbsFitterInfo(), jacobian_(TMatrixD(5, 5)), noise_(TMatrixDSym(5)) {
    reset();      
  }
  
  GblFitterInfo2::GblFitterInfo2(const TrackPoint* trackPoint, const AbsTrackRep* rep) : AbsFitterInfo(trackPoint, rep), jacobian_(TMatrixD(5, 5)), noise_(TMatrixDSym(5)) {
    reset();
    
    // Retrieve jacobian from rep by default (may not be used ... at 1st point)
    unsigned int dim = rep->getDim();
    if (dim != 5)
      throw new genfit::Exception("GblFitterInfo2: Representation state is not 5D", __LINE__, __FILE__);
    rep->getForwardJacobianAndNoise(jacobian_, noise_);      
    
  }
 
  GblFitterInfo2::GblFitterInfo2(const TrackPoint* trackPoint, const AbsTrackRep* rep, StateOnPlane& referenceState) : AbsFitterInfo(trackPoint, rep), jacobian_(TMatrixD(5, 5)) {
    reset();
    
    // Retrieve jacobian from rep by default (may not be used ... at 1st point)
    unsigned int dim = rep->getDim();
    if (dim != 5)
      throw new genfit::Exception("GblFitterInfo2: Representation state is not 5D", __LINE__, __FILE__);
    rep->getForwardJacobianAndNoise(jacobian_, noise_);
        
    setReferenceState(referenceState);
  }
  
  void GblFitterInfo2::reset(unsigned int repDim) {
    label_ = 0;
    measType_ = _noMeas;
    refPrediction_.ResizeTo(repDim);
    refPrediction_.Zero();
    
    TMatrixDSym zero(repDim);
    zero.Zero();
    
    fwdStateCorrection_.ResizeTo(repDim);
    fwdCov_.ResizeTo(zero);
    bwdStateCorrection_.ResizeTo(repDim);
    bwdCov_.ResizeTo(zero);
    diffFwdBwdPrediction_.ResizeTo(repDim);
    diffFwdBwdPrediction_.Zero();
    fwdPrediction_.ResizeTo(repDim);
    bwdPrediction_.ResizeTo(repDim);
    fwdCov_.SetMatrixArray(zero.GetMatrixArray());
    bwdCov_.SetMatrixArray(zero.GetMatrixArray());  
    
    measUsed_.clear();
    measDownWeight_.clear();
    measurement_.clear();
    measCov_.clear();
    hMatrix_.ResizeTo(HMatrixUV().getMatrix());
    hMatrix_ = HMatrixUV().getMatrix();
  }
  
  void GblFitterInfo2::setReferenceState(StateOnPlane& referenceState) {
    updateMeasurementAndPlane(referenceState);
    
    refPrediction_ = referenceState.getState();
    fwdPrediction_ = referenceState.getState();
    bwdPrediction_ = referenceState.getState();
    
    //TODO: reset even if already fitted?      
    fittedStateBwd_.reset(new MeasuredStateOnPlane(referenceState.getState(), trackPoint_->getTrack()->getCovSeed(), sharedPlane_, rep_, referenceState.getAuxInfo()));
    fittedStateFwd_.reset(new MeasuredStateOnPlane(referenceState.getState(), trackPoint_->getTrack()->getCovSeed(), sharedPlane_, rep_, referenceState.getAuxInfo()));
    
  }
  
  
  gbl::GblPoint GblFitterInfo2::constructGblPoint(unsigned int index, bool allowAmbiguities) {
    // All meas/scat info is added from genfit data again (to cope with possible RecoHit update)
    
    // Create GBL point with current jacobian
    gbl::GblPoint thePoint(jacobian_);    
    
    // Scatterer (all except first point)
    if (index > 0) {
      TMatrixDSym scatCov(4);
      // assuming no correlation between energy and scattering variance
      noise_.GetSub(1,4,scatCov, "S");
      thePoint.addThickScatterer(getKinksAndSteps(), scatCov.Invert());
    }
    
    // (implemented type of) Measurement
    if (measType_ != _noMeas) {
      unsigned int numMeasAllowed = allowAmbiguities ? measUsed_.size() : 1;
      for (unsigned int iMeas = 0; iMeas < numMeasAllowed; iMeas++) {
        TVectorD aRes(measurement_.at(iMeas) - hMatrix_ * fwdPrediction_);
        TMatrixDSym aCov(measCov_.at(iMeas));
        if (measType_ < _uvMeas) {
          // (1D) U or V measurement
          TVectorD aResidual(2);
          aResidual(measType_) = aRes(0);
          aResidual(1-measType_) = 0.;
          TVectorD aPrecision(2);
          aPrecision(measType_) = 1. / aCov(0,0);
          aPrecision(1-measType_) = 0.; // disable other measurement
          thePoint.addMeasurement(aResidual, aPrecision);
        } else {
          // (at least 2D) U,V measurement, arbitrary precision matrix
          thePoint.addMeasurement(aRes, aCov.Invert());
        }
        // init downweight
        measDownWeight_.at(measUsed_.at(iMeas)) = 1.;

        // Derivatives
        ICalibrationParametersDerivatives* globals = nullptr;
        if ((globals = dynamic_cast<ICalibrationParametersDerivatives*>(trackPoint_->getRawMeasurement(0)) )) {
          StateOnPlane sop(getFittedState().getState(), sharedPlane_, rep_, getFittedState(true).getAuxInfo());
          std::pair<std::vector<int>, TMatrixD> labelsAndMatrix = globals->globalDerivatives(&sop);
          std::vector<int> labels = labelsAndMatrix.first;
          TMatrixD derivs = labelsAndMatrix.second;

          if (derivs.GetNcols() > 0 && !labels.empty() && (unsigned int)derivs.GetNcols() == labels.size()) {
            thePoint.addGlobals(labels, derivs);
          }
          TMatrixD locals = globals->localDerivatives(&sop);
          if (locals.GetNcols() > 0) {
            thePoint.addLocals(locals);
            GblFitStatus* gblfs = dynamic_cast<GblFitStatus*>(trackPoint_->getTrack()->getFitStatus(rep_));
            if (gblfs) {
              if (gblfs->getMaxLocalFitParams() < locals.GetNcols())
                gblfs->setMaxLocalFitParams(locals.GetNcols());
            }
          }
        }
      }
    }
    
    return thePoint;
    
  }
  
  void GblFitterInfo2::updateMeasurementAndPlane(const StateOnPlane & sop) {
    if (!trackPoint_)
      return;
    if (!trackPoint_->hasRawMeasurements()) {
      //no measurement, plane is supposed to come from state (is perpendicular)
      setPlane(sop.getPlane());
      return;
    }
    std::vector<MeasurementOnPlane*> allMeas = trackPoint_->getRawMeasurement(0)->constructMeasurementsOnPlane(sop);

    // initial call ?
    bool initial = measDownWeight_.empty();
    if (initial) {
      // reset GBL downweights
      for (unsigned int i = 0; i < allMeas.size(); i++) {measDownWeight_.push_back(0.);}
      // get most probable measurement
      unsigned int imop = 0;
      double maxWeight = allMeas.at(0)->getWeight();
      for (unsigned int i = 1; i < allMeas.size(); i++)
        if (allMeas.at(i)->getWeight() > maxWeight) {
          imop = i;
          maxWeight = allMeas.at(i)->getWeight();
        }
      // use that
      measUsed_.push_back(imop);
    }
      
    unsigned int imop = measUsed_.at(0);
    hMatrix_.ResizeTo(allMeas.at(imop)->getHMatrix()->getMatrix());
    hMatrix_ = allMeas.at(imop)->getHMatrix()->getMatrix();
    // (implemented) measurement type(s)
    if (HMatrixU().getMatrix() == hMatrix_) { measType_ = _uMeas; }
    else if (HMatrixV().getMatrix() == hMatrix_) { measType_ = _vMeas; }
    else if (HMatrixUV().getMatrix() == hMatrix_) { measType_ = _uvMeas; }
    else {
      //no implemented measurement type, plane is supposed to come from state (is perpendicular)
      #ifdef DEBUG    
      std::cout << " no implemented meas type " << std::endl;
      #endif
      setPlane(sop.getPlane());
      return;
    }
    
    // check for ambiguities    
    if (initial) {  
      // ambiguity (CDC) ?
      if ((allMeas.size() == 2) && (measType_ == _uMeas)) { 
        // use that too
        measUsed_.push_back(1-imop);
      } 
    }
    
    // measurements to be used
    for (unsigned int i = 0; i < measUsed_.size(); i++) {
      unsigned int index = measUsed_.at(i);
      if (initial) {
        measurement_.push_back(allMeas.at(index)->getState());
        measCov_.push_back(allMeas.at(index)->getCov());
      } else {
        measurement_.at(i) = allMeas.at(index)->getState();
        measCov_.at(i) = allMeas.at(index)->getCov();
      }
    }  
    
    setPlane(allMeas.at(imop)->getPlane());
    
    for (unsigned int imeas = 0; imeas < allMeas.size(); imeas++)
      delete allMeas[imeas];
    allMeas.clear();
    
  }
  
  void GblFitterInfo2::resolveAmbiguities(gbl::GblTrajectory& traj, gbl::GblPoint& point, unsigned int &nResolved, unsigned int &nSwapped) {  
    if (!traj.isValid())
      return;
     
    // Implemented: one or two u or v measurement (-> max 2 meas. -> size of mResiduals etc)
    
    // (1 or 2) 1D u or v measurement
    if (measType_ < _uvMeas) {
      unsigned int numMRes = 0;        
      TVectorD mResiduals(2), mMeasErrors(2), mResErrors(2), mDownWeights(2);    
      if (0 != traj.getMeasResults(label_, numMRes, mResiduals, mMeasErrors, mResErrors, mDownWeights))
        throw genfit::Exception(" NO measurement results ", __LINE__, __FILE__);
      
      // resolve side
      if (numMRes == 2) {
        if (fabs(mResiduals[0]) <= fabs(mResiduals[1]))
          (point.getMeasBegin() + 1)->setEnabled(false); // disable measurement 1
        else {
          point.getMeasBegin()->setEnabled(false); // disable measurement 0
          // swap measurements (to restore order)
          unsigned int imop = measUsed_.at(0);
          measUsed_.at(0) = 1-imop;
          measUsed_.at(1) = imop;
          nSwapped++;
          #ifdef DEBUG    
          std::cout << " swap " << label_ << " " << mResiduals[0] << " " << mResiduals[1] << std::endl;
          #endif    
        }
        nResolved++;
      }
    }
    return;   
  }
  
  void GblFitterInfo2::updateDownweights(gbl::GblTrajectory& traj) {  
    if (!traj.isValid())
      return;
     
    // Implemented: single uv or one or two u or v measurement (-> max 2 meas. -> size of mResiduals etc)
    unsigned int numMRes = 0;        
    TVectorD mResiduals(2), mMeasErrors(2), mResErrors(2), mDownWeights(2);    
    if (0 != traj.getMeasResults(label_, numMRes, mResiduals, mMeasErrors, mResErrors, mDownWeights))
      throw genfit::Exception(" NO measurement results ", __LINE__, __FILE__);
    
    // (1 or 2) 1D u or v measurement
    if (measType_ < _uvMeas) {
      for (unsigned int i = 0; i < numMRes; ++i) {
        measDownWeight_.at(measUsed_.at(i)) = mDownWeights(i);
      }
      // reset weights for discarded measurements
      for (unsigned int i = numMRes; i < measUsed_.size(); ++i) {
        measDownWeight_.at(measUsed_.at(i)) = 0.;
      }
    }
    // (single) 2D uv measurement (2 downweights, use product)
    if (measType_ == _uvMeas && numMRes == 2) {
      measDownWeight_.at(measUsed_.at(0)) = mDownWeights(0)*mDownWeights(1);
    }
  }
    
  void GblFitterInfo2::updateFitResults(gbl::GblTrajectory& traj) {  
    if (!traj.isValid())
      return;
      
    unsigned int label = label_;  
    
    // Now update the the fitter info
    //
    //-------------------------------------------------
    // Backward/forward prediction (residual) (differs at scatterers) AFTER GBL fit
    
    // getResults respects size of ROOT matrices while converting from (internal) EIGEN3 matrices
    // forward prediction
    if (0 != traj.getResults(label, fwdStateCorrection_, fwdCov_))
      throw genfit::Exception(" NO forward results ", __LINE__, __FILE__);        
    
    // backward prediction
    if (0 != traj.getResults(-1 * label, bwdStateCorrection_, bwdCov_))
      throw genfit::Exception(" NO backward results ", __LINE__, __FILE__);        
    
    diffFwdBwdPrediction_ += fwdStateCorrection_ - bwdStateCorrection_;
    fwdPrediction_ += fwdStateCorrection_; // This is the update!
    bwdPrediction_ += bwdStateCorrection_; // This is the update!
    
    fittedStateFwd_.reset( new MeasuredStateOnPlane(fwdPrediction_, fwdCov_, sharedPlane_, rep_, getFittedState(true).getAuxInfo()) );         
    fittedStateBwd_.reset( new MeasuredStateOnPlane(bwdPrediction_, bwdCov_, sharedPlane_, rep_, getFittedState(true).getAuxInfo()) );
    
    //----------------------------------------------

    #ifdef DEBUG    
    // Residuals (in plane)
    // Implemented: single uv or one or two u or v measurement (-> max 2 meas. -> size of mResiduals etc)
    unsigned int numMRes = 0;    
    TVectorD mResiduals(2), mMeasErrors(2), mResErrors(2), mDownWeights(2);    
    if (0 != traj.getMeasResults(label, numMRes, mResiduals, mMeasErrors, mResErrors, mDownWeights))
      throw genfit::Exception(" NO measurement results ", __LINE__, __FILE__); 
    
    // Kinks and steps
    unsigned int numKRes = 0;
    TVectorD kResiduals(4), kMeasErrors(4), kResErrors(4), kDownWeights(4);    
    if (0 != traj.getScatResults(label, numKRes, kResiduals, kMeasErrors, kResErrors, kDownWeights))
      throw genfit::Exception(" NO scattering results ", __LINE__, __FILE__);        

    std::cout << std::endl;
    std::cout << " GBL residuals at label " << label << std::endl; 
    for (unsigned int i = 0; i < numMRes; ++i) {
      std::cout << " meas " << std::setw(15) << mResiduals(i) << " " << std::setw(15) << mMeasErrors(i) << " " 
      << std::setw(15) << mResErrors(i) << " " << std::setw(15) << mDownWeights(i) << std::endl;
    }     
    for (unsigned int i = 0; i < numKRes; ++i) {
      std::cout << " scat " << std::setw(15) << kResiduals(i) << " " << std::setw(15) << kMeasErrors(i) << " " 
      << std::setw(15) << kResErrors(i) << " " << std::setw(15) << kDownWeights(i) << std::endl;
    } 
    std::cout << std::endl;
    #endif
    
    //-------------------------------------------------
  }
  
  void GblFitterInfo2::recalculateJacobian(GblFitterInfo2* prevFitterInfo)
  {
    // Invalidates errors and corrections from last iteration
    // (will be defined in different plane). But fitted state and residual is ok.
        
    if (!prevFitterInfo) {
      jacobian_.UnitMatrix();
      noise_.Zero();
      
      prevFitterInfo = this;    
      return;
    }    
    
    //TODO
    //updateMeasurementAndPlane(StateOnPlane(fwdPrediction_, sharedPlane_, rep_));
    //
    
    // Take forward state from previous fitter info,
    // its (maybe updated) plane
    // and our rep
    StateOnPlane prevState(prevFitterInfo->getFittedState(true).getState(), prevFitterInfo->getPlane(), rep_, getFittedState(true).getAuxInfo());
    
    if (hasMeasurements()) {
      SharedPlanePtr newPlane = trackPoint_->getRawMeasurement(0)->constructPlane(prevState);
       rep_->extrapolateToPlane(prevState, newPlane, false, true);
    } else {
      rep_->extrapolateToPlane(prevState, sharedPlane_, false, true);      
    }
    
    rep_->getForwardJacobianAndNoise(jacobian_, noise_);
    // Now update meas data
    updateMeasurementAndPlane(prevState);
    
    //
    // Extrapolate predictions to new plane
    //
    bwdPrediction_ = prevState.getState();
    fwdPrediction_ = bwdPrediction_ + diffFwdBwdPrediction_;
    fittedStateFwd_.reset(new MeasuredStateOnPlane(fwdPrediction_, fwdCov_, sharedPlane_, rep_, getFittedState(true).getAuxInfo()));
    fittedStateBwd_.reset(new MeasuredStateOnPlane(bwdPrediction_, bwdCov_, sharedPlane_, rep_, getFittedState(true).getAuxInfo()));
    //
  }

  
  const MeasuredStateOnPlane& GblFitterInfo2::getFittedState(bool afterKink) const {      
    // ALways biased from GBL (global fit!)
    
    if (!fittedStateFwd_ || !fittedStateBwd_) { 
      //NOTE: This should be already set (from reference)! The auxInfo is being book-kept by it. If reference is not set, default auxInfo is used
      fittedStateFwd_.reset(new MeasuredStateOnPlane(fwdPrediction_, fwdCov_, sharedPlane_, rep_));        
      fittedStateBwd_.reset(new MeasuredStateOnPlane(bwdPrediction_, bwdCov_, sharedPlane_, rep_));
    }
    
    if (afterKink) {
      return *fittedStateFwd_;        
    }
    else {   
      return *fittedStateBwd_;        
    }
    
  }
  
  MeasurementOnPlane GblFitterInfo2::getResidual(unsigned int iMeas, bool, bool) const { 
    // only measurement errors implemented:
    //TODO: 1D covariance + more dimensions of residuals
    
    if (hasMeasurements()){
      // this means only for reference state before gbl fit, this way will be used
      TVectorD res( measurement_.at(iMeas) - hMatrix_ * fwdPrediction_ );      
      return MeasurementOnPlane(res, measCov_.at(iMeas), sharedPlane_, rep_, trackPoint_->getRawMeasurement(0)->constructHMatrix(getRep()));    
    }
    TMatrixDSym zeroCov(2);
    zeroCov.Zero();
    TVectorD zeroRes(2);
    zeroRes.Zero();
    // Else return 0's or whatever
    //TODO: or throw?
    return MeasurementOnPlane(zeroRes, zeroCov, sharedPlane_, rep_, new HMatrixUV());
    
  } // calculate residual
  
  TVectorD GblFitterInfo2::getKinksAndSteps() const {
    // return kinks and steps residuals = measurement(=0.) - prediction
    TVectorD kinksAndSteps(4);
    for (unsigned int i = 0; i < 4; ++i) { kinksAndSteps(i)= -diffFwdBwdPrediction_(i+1); }
    return kinksAndSteps;
  }
  
  MeasurementOnPlane GblFitterInfo2::getMeasurement(unsigned int iMeas) const{
    return MeasurementOnPlane(measurement_.at(iMeas), measCov_.at(iMeas), sharedPlane_, rep_, hasMeasurements() ? trackPoint_->getRawMeasurement(0)->constructHMatrix(rep_) : new HMatrixUV() );
  }
  
  bool GblFitterInfo2::checkConsistency(const genfit::PruneFlags*) const {
    //TODO
    return true;      
  }
  
  GblFitterInfo2* GblFitterInfo2::clone() const {
    
    GblFitterInfo2* retVal = new GblFitterInfo2(this->getTrackPoint(), this->getRep());
    
    retVal->setPlane(sharedPlane_);
    
    retVal->measType_ = measType_;
    retVal->label_ = label_;
    retVal->jacobian_ = jacobian_;
    retVal->noise_ = noise_;
    retVal->bwdStateCorrection_ = bwdStateCorrection_;
    retVal->fwdStateCorrection_ = fwdStateCorrection_;
    retVal->bwdCov_ = bwdCov_;
    retVal->fwdCov_ = fwdCov_;
    retVal->diffFwdBwdPrediction_ = diffFwdBwdPrediction_;
    retVal->fwdPrediction_ = fwdPrediction_;
    retVal->bwdPrediction_ = bwdPrediction_;
    retVal->refPrediction_ = refPrediction_;
    retVal->measUsed_ = measUsed_;
    retVal->measDownWeight_ = measDownWeight_;
    retVal->measurement_ = measurement_;
    retVal->measCov_ = measCov_;
    retVal->hMatrix_.ResizeTo(hMatrix_);
    retVal->hMatrix_ = hMatrix_;
    
    return retVal;
  }
  
  void GblFitterInfo2::Print(const Option_t*) const {
    //TODO
    std::cout << "=============================================================================================" << std::endl;
    std::cout << " >>>  GblFitterInfo2 " << std::endl;
    std::cout << "      ************** " << std::endl;
    
    std::cout << "                      rep: " << rep_ << ", trackpoint: " << trackPoint_ << ", plane: " << sharedPlane_.get() << std::endl;
    sharedPlane_->Print();
    std::cout << std::endl;
    
    std::cout << "===========================================================================================" << std::endl;
    std::cout << "     |          PREDICTIONS            |   REFERENCE    | Corrections from last iteration |" << std::endl;
    std::cout << "     | (+)prediction  | (-)prediction  |     state      | (+)correction |  (-) correction |" << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------" << std::endl;
    
    for (int i = 0; i <5; i++) {
      std::cout << std::left;
      std::cout << " ";
      if (i==0)
        std::cout << "q/p";
      if (i==1)
        std::cout << "u' ";
      if (i==2)
        std::cout << "v' ";
      if (i==3)
        std::cout << "u  ";
      if (i==4)
        std::cout << "v  ";
      std::cout << " | " 
      << std::setw(12) << fwdPrediction_(i)              << "   | "
      << std::setw(12) << bwdPrediction_(i)              << "   | "
      << std::setw(12) << refPrediction_(i)              << "   | "
      << std::setw(12) << fwdStateCorrection_(i)           << "   | "
      << std::setw(12) << bwdStateCorrection_(i)           << "   | " << std::endl;
    }
    std::cout << "===========================================================================================" << std::endl;    
    
    TVectorD residual(getKinksAndSteps());
    std::cout << "     | Scat. res.  u' | Scat. res.  v' | Scat. res.  u  | Scat. res.  v  |" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << "     | " 
    << std::setw(12) << residual(0) << "   | " 
    << std::setw(12) << residual(1) << "   | "
    << std::setw(12) << residual(2) << "   | "
    << std::setw(12) << residual(3) << "   | "
    << std::endl;
    
    std::cout << "=============================================================================================" << std::endl;
    std::cout << "Measurement type (u,v,uv,none): " << measType_ << std::endl;
    for (unsigned int iMeas = 0; iMeas < measurement_.size(); iMeas++) {
      std::cout << "Measurement : "; measurement_.at(iMeas).Print();
      std::cout << "H Matrix : "; hMatrix_.Print();
      std::cout << "Measurement covariance : "; measCov_.at(iMeas).Print();
      std::cout << "Measurement residuals  : "; getResidual(iMeas).getState().Print();
      std::cout << "Measurement down-weight: " << measDownWeight_.at(iMeas) << std::endl;
    }
    std::cout << "Jacobian: "; jacobian_.Print();
    std::cout << "Noise   : "; noise_.Print();
    std::cout << "Backward covariance: "; bwdCov_.Print();
    std::cout << "Forward covariance : "; fwdCov_.Print();

    std::cout << "=============================================================================================" << std::endl;    
    
  }

  
} // end of namespace genfit
