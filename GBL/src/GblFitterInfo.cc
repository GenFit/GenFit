#include "GblFitterInfo.h"
#include "MeasurementOnPlane.h"
#include "HMatrixUV.h"
#include "HMatrixV.h"
#include "HMatrixU.h"

//! Namespace for the general broken lines package
namespace genfit {

  GblFitterInfo::GblFitterInfo() : AbsFitterInfo(), jacobian_(TMatrixD(5, 5)) {
    reset();      
  }
  
  GblFitterInfo::GblFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep) : AbsFitterInfo(trackPoint, rep), jacobian_(TMatrixD(5, 5)) {
    reset();
    
    // Retrieve jacobian from rep by default (may not be used ... at 1st point)
    unsigned int dim = rep->getDim();
    if (dim != 5)
      throw new genfit::Exception("GblFitterInfo: Representation state is not 5D", __LINE__, __FILE__);
    TMatrixDSym noise(dim, dim);
    TVectorD dState(dim);
    rep->getForwardJacobianAndNoise(jacobian_, noise, dState);      
    
  }
 
  GblFitterInfo::GblFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep, StateOnPlane& referenceState) : AbsFitterInfo(trackPoint, rep), jacobian_(TMatrixD(5, 5)) {
    reset();
    
    // Retrieve jacobian from rep by default (may not be used ... at 1st point)
    unsigned int dim = rep->getDim();
    if (dim != 5)
      throw new genfit::Exception("GblFitterInfo: Representation state is not 5D", __LINE__, __FILE__);
    TMatrixDSym noise(dim, dim);
    TVectorD dState(dim);
    rep->getForwardJacobianAndNoise(jacobian_, noise, dState);
        
    setReferenceState(referenceState);
  }
  
  void GblFitterInfo::reset(unsigned int measurementDim, unsigned int repDim) {
    refPrediction_.ResizeTo(repDim);
    refPrediction_.Zero();
    
    measResiduals_.ResizeTo(measurementDim);
    measResiduals_.Zero();
    kinkResiduals_.ResizeTo(measurementDim);
    kinkResiduals_.Zero();
    measResidualErrors_.ResizeTo(measurementDim);
    measResidualErrors_.Zero();
    kinkResidualErrors_.ResizeTo(measurementDim);
    kinkResidualErrors_.Zero();
    measDownWeights_.ResizeTo(measurementDim);
    measDownWeights_.Zero();
    kinkDownWeights_.ResizeTo(measurementDim);
    kinkDownWeights_.Zero();
    
    TMatrixDSym zero(repDim);
    zero.Zero();
    
    fwdStateCorrection_.ResizeTo(repDim);
    fwdCov_.ResizeTo(zero);
    bwdStateCorrection_.ResizeTo(repDim);
    bwdCov_.ResizeTo(zero);
    fwdPrediction_.ResizeTo(repDim);
    bwdPrediction_.ResizeTo(repDim);
    fwdCov_.SetMatrixArray(zero.GetMatrixArray());
    bwdCov_.SetMatrixArray(zero.GetMatrixArray());  
    
    measurement_.ResizeTo(2);
    measurement_.Zero();
    measCov_.ResizeTo(TMatrixDSym(2));
    measCov_.Zero();
    hMatrix_.ResizeTo(HMatrixUV().getMatrix());
    hMatrix_ = HMatrixUV().getMatrix();
  }
  
  void GblFitterInfo::setReferenceState(StateOnPlane& referenceState) {
    updateMeasurementAndPlane(referenceState);
    
    refPrediction_ = referenceState.getState();
    fwdPrediction_ = referenceState.getState();
    bwdPrediction_ = referenceState.getState();
    
    //TODO: reset even if already fitted?      
    fittedStateBwd_.reset(new MeasuredStateOnPlane(referenceState, trackPoint_->getTrack()->getCovSeed()));
    fittedStateFwd_.reset(new MeasuredStateOnPlane(referenceState, trackPoint_->getTrack()->getCovSeed()));
    
  }
  
  void GblFitterInfo::setJacobian(TMatrixD jacobian) {
    jacobian_.ResizeTo(jacobian);
    jacobian_ = jacobian; 
  }
  
  TMatrixDSym GblFitterInfo::getCovariance(double variance, TVector3 trackDirection, SharedPlanePtr measurementPlane) const {
    
    double c1 = trackDirection.Dot(measurementPlane->getU());
    double c2 = trackDirection.Dot(measurementPlane->getV());
    
    TMatrixDSym scatCov(2);
    scatCov(0, 0) = 1. - c2 * c2;
    scatCov(1, 1) = 1. - c1 * c1;
    scatCov(0, 1) = c1 * c2;
    scatCov(1, 0) = c1 * c2;
    scatCov *= variance * variance / (1. - c1 * c1 - c2 * c2) / (1. - c1 * c1 - c2 * c2) ;
    
    return scatCov;
  }
  
  gbl::GblPoint GblFitterInfo::constructGblPoint() {
    // All meas/scat info is added from genfit data again (to cope with possible RecoHit update)
    
    // Create GBL point with current jacobian
    gbl::GblPoint thePoint(jacobian_);
    
    //NOTE: 3rd update and update anytime GblPoint is requested
    // mostly likely will update with reference as on 2nd update
    StateOnPlane sop(getFittedState().getState(), sharedPlane_, rep_);
    
    
    // Scatterer
    if (trackPoint_->hasThinScatterer()) {
      if (!hasMeasurements()) {
        //double scatVariance = trackPoint_->getMaterialInfo()->getMaterial().getDensity();
        //TVectorD kinkPrec(2);
        //kinkPrec(0) = 1./scatVariance/scatVariance; kinkPrec(1) = 1./scatVariance/scatVariance;
        //thePoint.addScatterer(getKinks(), kinkPrec);
        //TODO: if state at scatterer is updated, the direction of track might not be perpendicular anymore
        // plane does not change at pure scatterer
        TMatrixDSym kinkCov = getCovariance(trackPoint_->getMaterialInfo()->getMaterial().density, sop.getDir(), sop.getPlane());
        thePoint.addScatterer(getKinks(), kinkCov.Invert());
      } else {
        TMatrixDSym kinkCov = getCovariance(trackPoint_->getMaterialInfo()->getMaterial().density, sop.getDir(), sop.getPlane());
        thePoint.addScatterer(getKinks(), kinkCov.Invert());
      }      
    }
    
    
    
    // Measurement
    if (hasMeasurements()) {
      MeasuredStateOnPlane measurement = getResidual(0, true, true);    
      TVectorD aResiduals(measurement.getState());
      TMatrixDSym aPrecision(measurement.getCov().Invert()); 
      if (HMatrixU().getMatrix() == hMatrix_) {
        double res = aResiduals(0);
        double prec = aPrecision(0, 0);
        aResiduals.ResizeTo(2);
        aPrecision.ResizeTo(TMatrixDSym(2));
        aResiduals.Zero();
        aResiduals(0) = res;
        aPrecision.Zero();
        aPrecision(0, 0) = prec;
      }      
      if (HMatrixV().getMatrix() == hMatrix_) {
        double res = aResiduals(0);
        double prec = aPrecision(0, 0);
        aResiduals.ResizeTo(2);
        aPrecision.ResizeTo(TMatrixDSym(2));
        aResiduals.Zero();
        aResiduals(1) = res;
        aPrecision.Zero();
        aPrecision(1, 1) = prec;
      }
      // always 2D, U/V set by precision (=0 to disable)
      thePoint.addMeasurement(aResiduals, aPrecision);
    }
    
    // Derivatives      
    ICalibrationParametersDerivatives* globals = nullptr;
    if (hasMeasurements() && (globals = dynamic_cast<ICalibrationParametersDerivatives*>(trackPoint_->getRawMeasurement(0)) )) {    
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
    
    return thePoint;
    
  }
  
  void GblFitterInfo::updateMeasurementAndPlane(const StateOnPlane & sop) {
    if (!trackPoint_)
      return;
    if (!trackPoint_->hasRawMeasurements()) {
      //no measurement, plane is supposed to come from state (is perpendicular)
      setPlane(sop.getPlane());
      return;
    }
    std::vector<MeasurementOnPlane*> allMeas = trackPoint_->getRawMeasurement(0)->constructMeasurementsOnPlane(sop);
    
    /*
    double normMin(9.99E99);
    unsigned int imop(0);
    const AbsHMatrix* H = allMeas[0]->getHMatrix();
    for (unsigned int i=0; i<allMeas.size(); ++i) {
      if (*(allMeas[i]->getHMatrix()) != *H){
        Exception e("GblFitterInfo::updateMeasurementAndPlane: Cannot compare measurements with different H-Matrices.", __LINE__,__FILE__);
        e.setFatal();
        throw e;
      }
      
      TVectorD res = allMeas[i]->getState() - H->Hv(sop.getState());
      double norm = sqrt(res.Norm2Sqr());
      if (norm < normMin) {
        normMin = norm;
        imop = i;
      }
    }
    */
    unsigned int imop = 0;
    double maxWeight = allMeas.at(0)->getWeight();
    for (unsigned int i = 0; i < allMeas.size(); i++)
      if (allMeas.at(i)->getWeight() > maxWeight)
        imop = i;

    measurement_.ResizeTo(allMeas.at(imop)->getState());
    measurement_ = allMeas.at(imop)->getState();
    measCov_.ResizeTo(allMeas.at(imop)->getCov());
    measCov_ = allMeas.at(imop)->getCov();
    hMatrix_.ResizeTo(allMeas.at(imop)->getHMatrix()->getMatrix());
    hMatrix_ = allMeas.at(imop)->getHMatrix()->getMatrix();
    
    setPlane(allMeas.at(imop)->getPlane());
    
    for (unsigned int imeas = 0; imeas < allMeas.size(); imeas++)
      delete allMeas[imeas];
    allMeas.clear();
    

  }
  
  void GblFitterInfo::updateFitResults(gbl::GblTrajectory& traj) {  
    if (!traj.isValid())
      return;
    
    // Deduce our position in the trajectory
    //----------------------------------------------
    unsigned int label = 0;
    // Loop over track and find index of this fitter info
    genfit::Track* trk = trackPoint_->getTrack();
    for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {
      // Indexing of GblFitter info in track (each gives one point to trajectory at trackpoint position)
      if (dynamic_cast<GblFitterInfo*>( trk->getPoint(ip)->getFitterInfo(rep_) )) {
        // First fitter info has label 1 (for point in GBL traj)
        label++;
        // We found itself = we have our label
        if (trk->getPoint(ip)->getFitterInfo(rep_) == this)
          break;
      }
    }
    if (label == 0)
      throw genfit::Exception("GblFitterInfo: fitter info did not found itself in track to update", __LINE__, __FILE__);      
    if (label > traj.getNumPoints())
      throw genfit::Exception("GblFitterInfo: Deduced point label not valid", __LINE__, __FILE__);    
    //----------------------------------------------
    
    // Residuals (in plane)
    unsigned int numMRes = 2;    
    TVectorD mResiduals(2), mMeasErrors(2), mResErrors(2), mDownWeights(2);    
    if (0 != traj.getMeasResults(label, numMRes, mResiduals, mMeasErrors, mResErrors, mDownWeights))
      throw genfit::Exception(" NO measurement results ", __LINE__, __FILE__);  
    
    // Kinks
    unsigned int numKRes = 2;
    TVectorD kResiduals(2), kMeasErrors(2), kResErrors(2), kDownWeights(2);    
    if (0 != traj.getScatResults(label, numKRes, kResiduals, kMeasErrors, kResErrors, kDownWeights))
      throw genfit::Exception(" NO scattering results ", __LINE__, __FILE__);        
    
    // Check for local derivatives    
    int nLocals = 0;
    GblFitStatus* gblfs = dynamic_cast<GblFitStatus*>(trackPoint_->getTrack()->getFitStatus(rep_));
    if (gblfs)
      nLocals = gblfs->getMaxLocalFitParams();
     
    // Predictions (different at scatterer)
    TVectorD bwdUpdate(5 + nLocals), fwdUpdate(5 + nLocals);
    TMatrixDSym bwdCov(5 + nLocals), fwdCov(5 + nLocals);
    
    // forward prediction
    if (0 != traj.getResults(label, fwdUpdate, fwdCov))
      throw genfit::Exception(" NO forward results ", __LINE__, __FILE__);        
    
    // backward prediction
    if (0 != traj.getResults(-1 * label, bwdUpdate, bwdCov))
      throw genfit::Exception(" NO backward results ", __LINE__, __FILE__);        
    
    if (nLocals > 0) {
      TVectorD _bwdUpdate(5 + nLocals), _fwdUpdate(5 + nLocals);
      TMatrixDSym _bwdCov(5 + nLocals), _fwdCov(5 + nLocals);
      _bwdUpdate = bwdUpdate;
      _fwdUpdate = fwdUpdate;
      _bwdCov = bwdCov;
      _fwdCov = fwdCov;
      bwdUpdate.ResizeTo(5);
      fwdUpdate.ResizeTo(5);
      bwdCov.ResizeTo(TMatrixDSym(5));
      fwdCov.ResizeTo(TMatrixDSym(5));
      for (int i = 0; i < 5; i++) {
        bwdUpdate(i) = _bwdUpdate(i);
        fwdUpdate(i) = _fwdUpdate(i);
        for (int j = 0; j < 5; j++) {
          bwdCov(i, j) = _bwdCov(i, j);
          fwdCov(i, j) = _fwdCov(i, j);
        }        
      }
    }
    // Now update the the fitter info
    //
    //-------------------------------------------------
    // Backward/forward prediction (residual) (differs at scatterers) AFTER GBL fit
    bwdStateCorrection_ = bwdUpdate;
    fwdStateCorrection_ = fwdUpdate;
    bwdCov_ = bwdCov;
    fwdCov_ = fwdCov;
    fwdPrediction_ += fwdStateCorrection_; // This is the update!
    bwdPrediction_ += bwdStateCorrection_; // This is the update!
    
    fittedStateFwd_.reset( new MeasuredStateOnPlane(fwdPrediction_, fwdCov_, sharedPlane_, rep_) );         
    fittedStateBwd_.reset( new MeasuredStateOnPlane(bwdPrediction_, bwdCov_, sharedPlane_, rep_) );
    
    // Set scattering/measurement residual data
    kinkResiduals_ = kResiduals;
    measResiduals_ = mResiduals;
    kinkResidualErrors_ = kResErrors;
    measResidualErrors_ = mResErrors;
    measDownWeights_ = mDownWeights;
    kinkDownWeights_ = kDownWeights;   
    
    //-------------------------------------------------
  }
  
  void GblFitterInfo::recalculateJacobian(GblFitterInfo* prevFitterInfo)
  {
    // Invalidates errors and corrections from last iteration
    // (will be defined in different plane). But fitted state and residual is ok.
        
    if (!prevFitterInfo) {
      jacobian_.UnitMatrix();
      
      //TODO: For 1st plane this has no sense
      //return;
      prevFitterInfo = this;
    }    
    
    //TODO
    SharedPlanePtr oldPlane(new genfit::DetPlane(*sharedPlane_));

    //updateMeasurementAndPlane(StateOnPlane(fwdPrediction_, sharedPlane_, rep_));
    //
    
    TMatrixDSym noise;
    TVectorD dstate;
    // Take forward state from previous fitter info,
    // its (maybe updated) plane
    // and our rep
    StateOnPlane prevState(prevFitterInfo->getFittedState(true).getState(), prevFitterInfo->getPlane(), rep_);
    
    if (hasMeasurements()) {
      SharedPlanePtr newPlane = trackPoint_->getRawMeasurement(0)->constructPlane(prevState);
      rep_->extrapolateToPlane(prevState, newPlane, false, true);
    } else {
      rep_->extrapolateToPlane(prevState, sharedPlane_, false, true);      
    }
    
    rep_->getForwardJacobianAndNoise(jacobian_, noise, dstate);
    // Now update meas data
    updateMeasurementAndPlane(prevState);
    
    //
    // Extrap predictions to new plane
    //
    StateOnPlane oldFwdState(fwdPrediction_, oldPlane, rep_);
    StateOnPlane oldBwdState(bwdPrediction_, oldPlane, rep_);
    rep_->extrapolateToPlane(oldFwdState, sharedPlane_);
    rep_->extrapolateToPlane(oldBwdState, sharedPlane_);
    fwdPrediction_ = oldFwdState.getState();
    bwdPrediction_ = oldBwdState.getState();
    fittedStateBwd_.reset();
    fittedStateFwd_.reset();
    //
  }

  
  const MeasuredStateOnPlane& GblFitterInfo::getFittedState(bool afterKink) const {      
    // ALways biased from GBL (global fit!)
    
    if (!fittedStateFwd_ || !fittedStateBwd_) { 
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
  
  MeasurementOnPlane GblFitterInfo::getResidual(unsigned int, bool, bool onlyMeasurementErrors) const { 
    // Return error of resduals (=0, none, for reference track, no errors known before fit (seed covariance might not be correct)
    TMatrixDSym localCovariance(2);
    localCovariance.Zero();
    // combined meas+fit errors:
    //TODO: 1D covariance + more dimensions of residuals
    
    localCovariance(0, 0) = measResidualErrors_(0) * measResidualErrors_(0);
    localCovariance(1, 1) = measResidualErrors_(1) * measResidualErrors_(1); 
    
    if (HMatrixU().getMatrix() == hMatrix_) {
      localCovariance.ResizeTo(TMatrixDSym(1));
      localCovariance(0, 0) = measResidualErrors_(0) * measResidualErrors_(0);
    }
    if (HMatrixV().getMatrix() == hMatrix_) {
      localCovariance.ResizeTo(TMatrixDSym(1));
      localCovariance(0, 0) = measResidualErrors_(1) * measResidualErrors_(1);
    }
    
    if (hasMeasurements()){
      // this means only for reference state before gbl fit, this way will be used
      TVectorD res( measurement_ - hMatrix_ * getFittedState(true).getState() );      
      
      if (onlyMeasurementErrors) {
        localCovariance.ResizeTo(measCov_);
        localCovariance = measCov_;
      }
      return MeasurementOnPlane(res, localCovariance, sharedPlane_, rep_, trackPoint_->getRawMeasurement(0)->constructHMatrix(getRep()));    
    }
    TVectorD zeroRes(2);
    zeroRes.Zero();
    // Else return 0's or whatever
    //TODO: or throw?
    return MeasurementOnPlane(zeroRes, localCovariance, sharedPlane_, rep_, new HMatrixUV());
    
  } // calculate residual, track and measurement errors are added if onlyMeasurementErrors is false
  
  MeasurementOnPlane GblFitterInfo::getKink() const { 
    TMatrixDSym localCovariance(2);
    localCovariance.Zero();
    localCovariance(0, 0) = kinkResidualErrors_(0) * kinkResidualErrors_(0);
    localCovariance(1, 1) = kinkResidualErrors_(1) * kinkResidualErrors_(1);      
    return MeasurementOnPlane(getKinks(), localCovariance, sharedPlane_, rep_, new genfit::HMatrixUV());
    
  } 
  
  TVectorD GblFitterInfo::getKinks() const {
    TVectorD kinks(2);
    kinks.Zero();
    
    TVectorD stateDiff(getFittedState(true).getState() - getFittedState(false).getState());
    kinks(0) = -stateDiff(1);
    kinks(1) = -stateDiff(2);
    
    return kinks;      
  }
  
  MeasurementOnPlane GblFitterInfo::getMeasurement() const{
    return MeasurementOnPlane(measurement_, measCov_, sharedPlane_, rep_, hasMeasurements() ? trackPoint_->getRawMeasurement(0)->constructHMatrix(rep_) : new HMatrixUV() );
  }
  
  bool GblFitterInfo::checkConsistency(const genfit::PruneFlags*) const {
    //TODO
    return true;      
  }
  
  GblFitterInfo* GblFitterInfo::clone() const {
    
    GblFitterInfo* retVal = new GblFitterInfo(this->getTrackPoint(), this->getRep());
    
    retVal->setPlane(sharedPlane_);
    
    retVal->jacobian_ = jacobian_;
    retVal->measResiduals_ = measResiduals_;
    retVal->measResidualErrors_ = measResidualErrors_;
    retVal->kinkResiduals_ = kinkResiduals_;
    retVal->kinkResidualErrors_ = kinkResidualErrors_;
    retVal->measDownWeights_ = measDownWeights_;
    retVal->kinkDownWeights_ = kinkDownWeights_;
    retVal->bwdStateCorrection_ = bwdStateCorrection_;
    retVal->fwdStateCorrection_ = fwdStateCorrection_;
    retVal->bwdCov_ = bwdCov_;
    retVal->fwdCov_ = fwdCov_;
    retVal->fwdPrediction_ = fwdPrediction_;
    retVal->bwdPrediction_ = bwdPrediction_;
    retVal->refPrediction_ = refPrediction_;
    retVal->measurement_.ResizeTo(measurement_);
    retVal->measCov_.ResizeTo(measCov_);
    retVal->hMatrix_.ResizeTo(hMatrix_);
    retVal->measurement_ = measurement_;
    retVal->measCov_ = measCov_;
    retVal->hMatrix_ = hMatrix_;
    
    return retVal;
  }
  
  void GblFitterInfo::Print(const Option_t*) const {
    //TODO
    std::cout << "=============================================================================================" << std::endl;
    std::cout << " >>>  GblFitterInfo " << std::endl;
    std::cout << "      ************* " << std::endl;
    
    std::cout << "                      rep: " << rep_ << ", trackpoint: " << trackPoint_ << ", plane: " << sharedPlane_.get() << std::endl;
    sharedPlane_->Print();
    std::cout << std::endl;
    
    std::cout << "=============================================================================================" << std::endl;
    std::cout << "     |          PREDICTIONS            |   REFERENCE    |  Corrections from last iteration  |" << std::endl;
    std::cout << "     | (+)prediction  | (-)prediction  |     state      |  (+)correction |  (-) correction  |" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------" << std::endl;
    
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
      << std::setw(12) << bwdStateCorrection_(i)           << std::endl;
    }
    std::cout << "=============================================================================================" << std::endl;    
    
    if (hasMeasurements()) {
      std::cout << "     | Meas. residual |     measurement - prediction    |   Down-weight  |   Fit+meas Err.  |" << std::endl;
      std::cout << "     |                |                                 |                |   -diagonaliz.   |" << std::endl;
      std::cout << "---------------------------------------------------------------------------------------------" << std::endl;
      
      TVectorD residual = getResidual().getState();
      if (residual.GetNoElements()<2) {
        residual.ResizeTo(2);
        residual.Zero();
        if (trackPoint_->getRawMeasurement(0)->constructHMatrix(getRep())->isEqual(genfit::HMatrixU()))
          residual(0) = getResidual().getState()(0);
        else
          residual(1) = getResidual().getState()(0);
        
      }
      std::cout << " u   | " 
      << std::setw(12) << measResiduals_(0)      << "   |  "      
      << std::setw(12) << residual(0)            << "                   |  "
      << std::setw(12) << measDownWeights_(0)    << "  |  "
      << std::setw(12) << measResidualErrors_(0) 
      << std::endl;
      
      std::cout << " v   | " 
      << std::setw(12) << measResiduals_(1)      << "   |  "      
      << std::setw(12) << residual(1)            << "                   |  "
      << std::setw(12) << measDownWeights_(1)    << "  |  "
      << std::setw(12) << measResidualErrors_(1) 
      << std::endl;
      
      std::cout << "---------------------------------------------------------------------------------------------" << std::endl;
    }
    
    std::cout << "     | Kink residual  |  Residual of slope difference   |   Down-weight  |   Fit Kink Err.  |" << std::endl;
    std::cout << "     | -diagonalized  |   - (   (+)pred - (-)pred   )   |                |   -diagonaliz.   |" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------" << std::endl;
    
    std::cout << " u'  | " 
    << std::setw(12) << kinkResiduals_(0)        << "   |  "            
    << std::setw(12) << getKinks()(0)            << "                   |  "
    << std::setw(12) << kinkDownWeights_(0)      << "  |  "
    << std::setw(12) << kinkResidualErrors_(0) 
    << std::endl;
    
    std::cout << " v'  | " 
    << std::setw(12) << kinkResiduals_(1)        << "   |  "            
    << std::setw(12) << getKinks()(1)            << "                   |  "
    << std::setw(12) << kinkDownWeights_(1)      << "  |  "
    << std::setw(12) << kinkResidualErrors_(1) 
    << std::endl;
    std::cout << "=============================================================================================" << std::endl;    
    std::cout << "Measurement: "; measurement_.Print();
    
    std::cout << "H Matrix: "; hMatrix_.Print();
    
    std::cout << "Measurement covariance: "; measCov_.Print();
    
    std::cout << "Jacobian: "; jacobian_.Print();
    std::cout << "Backward covariance: "; bwdCov_.Print();
    std::cout << "Forward covariance : "; fwdCov_.Print();
    
    std::cout << "=============================================================================================" << std::endl;    
    
  }

  
} // end of namespace genfit
