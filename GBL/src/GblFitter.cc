//-*-mode: C++; c-basic-offset: 2; -*-
/* Copyright 2013-2014
 *  Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *  This is an interface to General Broken Lines
 * 
 *  Version: 6 --------------------------------------------------------------
 *  - complete rewrite to GblFitter using GblFitterInfo and GblFitStatus
 *  - mathematics should be the same except for additional iterations 
 *    (not possible before)
 *  - track is populated with scatterers + additional points with 
 *    scatterers and no measurement (optional)
 *  - final track contains GblFitStatus and fitted states from GBL prediction
 *  - 1D/2D hits supported (pixel, single strip, combined strips(2D), wire)
 *  - At point: Only the very first raw measurement is used and from
 *    that, constructed MeasurementOnPlane with heighest weight
 *  Version: 5 --------------------------------------------------------------
 *  - several bug-fixes:
 *    - Scatterers at bad points
 *    - Jacobians at a point before they should be (code reorganized)
 *    - Change of sign of residuals
 *  Version: 4 --------------------------------------------------------------
 *    Fixed calculation of equvivalent scatterers (solution by C. Kleinwort)
 *    Now a scatterer is inserted at each measurement (except last) and
 *    between each two measurements.
 *    TrueHits/Clusters. Ghost (1D) hits ignored. With or
 *    without magnetic field.
 *  Version: 3 --------------------------------------------------------------
 *    This version now supports both TrueHits and Clusters for VXD.
 *    It can be used for arbitrary material distribution between
 *    measurements. Moments of scattering distribution are computed
 *    and translated into two equivalent thin GBL scatterers placed
 *    at computed positions between measurement points.
 *  Version: 2 ... never published -----------------------------------------
 *    Scatterer at each boundary (tooo many scatterers). TrueHits/Clusters.
 *    Without global der.&MP2 output.
 *  Version: 1 --------------------------------------------------------------
 *    Scatterers at measurement planes. TrueHits
 *  Version: 0 --------------------------------------------------------------
 *    Without scatterers. Genfit 1.
 *  -------------------------------------------------------------------------
 *
 *  This file is part of GENFIT.
 *
 *  GENFIT is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GENFIT is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "GblFitter.h"
#include "../include/GblFitStatus.h"
#include "GblFitterInfo.h"
#include "GblTrajectory.h"
#include "GblPoint.h"
#include "ICalibrationParametersDerivatives.h"

#include "Track.h"
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <string>
#include <list>
#include <FieldManager.h>
#include <HMatrixU.h>
#include <HMatrixV.h>
#include <Math/SMatrix.h>
#include <TMatrixD.h>
#include <TVectorDfwd.h>
#include <TMatrixT.h>
#include <TVector3.h>

//#define DEBUG

using namespace gbl;
using namespace std;
using namespace genfit;

/**
 * Destructor
 */
GblFitter::~GblFitter() {
  if (m_segmentController) {
    delete m_segmentController;
    m_segmentController = nullptr;
  }
}

void GblFitter::setTrackSegmentController(GblTrackSegmentController* controler)
{
  if (m_segmentController) {
    delete m_segmentController;
    m_segmentController = nullptr;
  }
  m_segmentController = controler;      
}

void GblFitter::processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits)
{
  cleanGblInfo(trk, rep);
  
  if (resortHits)
    sortHits(trk, rep);
  
  // This flag enables/disables fitting of q/p parameter in GBL
  // It is switched off automatically if no B-field at (0,0,0) is detected.
  bool fitQoverP = true;
  //TODO: Use clever way to determine zero B-field
  double Bfield = genfit::FieldManager::getInstance()->getFieldVal(TVector3(0., 0., 0.)).Mag();
  if (!(Bfield > 1.e-16))
    fitQoverP = false;
  // degrees of freedom after fit
  int Ndf = 0;
  // Chi2 after fit
  double Chi2 = 0.;
  //FIXME: d-w's not used so far...
  double lostWeight = 0.;  

  // Preparation of points (+add reference states) for GBL fit
  // -----------------------------------------------------------------
  genfit::GblFitStatus* gblfs = new genfit::GblFitStatus();
  trk->setFitStatus(gblfs, rep);
  gblfs->setCurvature(fitQoverP);
  //
  // Propagate reference seed, create scattering points, calc Jacobians
  // and store everything in fitter infos. (ready to collect points and fit)
  //
  //
  gblfs->setTrackLen(
  //
    constructGblInfo(trk, rep)
  //
  );
  //
  //
  gblfs->setIsFittedWithReferenceTrack(true);
  gblfs->setNumIterations(0); //default value, still valid, No GBL iteration
  if (m_externalIterations < 1)
    return;
  // -----------------------------------------------------------------
  

  unsigned int nFailed = 0;
  int fitRes = 0;
  std::vector<std::string> gblIterations;
  gblIterations.push_back(m_gblInternalIterations);

  // Iterations and updates of fitter infos and fit status
  // ------------------------------------------------------------------- 
  for (unsigned int iIter = 0; iIter < m_externalIterations; iIter++) {
    // GBL refit (1st of reference, then refit of GBL trajectory itself)
    int nscat = 0, nmeas = 0, ndummy = 0;
    std::vector<gbl::GblPoint> points = collectGblPoints(trk, rep);
    for(unsigned int ip = 0;ip<points.size(); ip++) {
      GblPoint & p = points.at(ip);
      if (p.hasScatterer())
        nscat++;
      if (p.hasMeasurement())
        nmeas++;
      if(!p.hasMeasurement()&&!p.hasScatterer())
        ndummy++;
    }
    gbl::GblTrajectory traj(points, gblfs->hasCurvature());
    
    fitRes = traj.fit(Chi2, Ndf, lostWeight, (iIter == m_externalIterations - 1) ? m_gblInternalIterations : "");
    
    // Update fit results in fitterinfos
    updateGblInfo(traj, trk, rep);
    
    // This repropagates to get new Jacobians,
    // if planes changed, predictions are extrapolated to new planes
    if (m_recalcJacobians > iIter) {
      GblFitterInfo* prevFitterInfo = 0;
      GblFitterInfo* currFitterInfo = 0;
      for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {
        if (trk->getPoint(ip)->hasFitterInfo(rep) && (currFitterInfo = dynamic_cast<GblFitterInfo*>(trk->getPoint(ip)->getFitterInfo(rep)))) {

          currFitterInfo->recalculateJacobian(prevFitterInfo);
          prevFitterInfo = currFitterInfo;
        }
      }
    }
        
    gblfs->setIsFitted(true);
    gblfs->setIsFitConvergedPartially(fitRes == 0);
    nFailed = trk->getNumPointsWithMeasurement() - nmeas;
    gblfs->setNFailedPoints(nFailed);
    gblfs->setIsFitConvergedFully(fitRes == 0 && nFailed == 0);
    gblfs->setNumIterations(iIter + 1);
    gblfs->setChi2(Chi2);    
    gblfs->setNdf(Ndf);
    gblfs->setCharge(trk->getFittedState().getCharge());
    
    #ifdef DEBUG
    int npoints_meas = trk->getNumPointsWithMeasurement();  
    int npoints_all = trk->getNumPoints();
        
    cout << "-------------------------------------------------------" << endl;
    cout << "               GBL processed genfit::Track            " << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << " # Track Points       :  " << npoints_all  << endl;
    cout << " # Meas. Points       :  " << npoints_meas << endl;
    cout << " # GBL points all     :  " << traj.getNumPoints();
    if (ndummy)
      cout << " (" << ndummy << " dummy) ";
    cout << endl;
    cout << " # GBL points meas    :  " << nmeas << endl;
    cout << " # GBL points scat    :  " << nscat << endl;    
    cout << "-------------- GBL Fit Results ----------- Iteration  " << iIter+1 << " " << ((iIter < gblIterations.size()) ? gblIterations[iIter] : "") << endl;
    cout << " Fit q/p parameter    :  " << (gblfs->hasCurvature() ? ("True") : ("False")) << endl;
    cout << " Valid trajectory     :  " << ((traj.isValid()) ? ("True") : ("False")) << endl;
    cout << " Fit result           :  " << fitRes << "    (0 for success)" << endl;
    cout << " GBL track NDF        :  " << Ndf << "    (-1 for failure)" << endl;
    cout << " GBL track Chi2       :  " << Chi2 << endl;
    cout << " GBL track P-value    :  " << TMath::Prob(Chi2, Ndf) << endl;
    cout << "-------------------------------------------------------" << endl;
    #endif
    
  }  
  // -------------------------------------------------------------------

}

void GblFitter::cleanGblInfo(Track* trk, const AbsTrackRep* rep) const {
  
  for (int ip = trk->getNumPoints() - 1; ip >=0; ip--) {
    trk->getPoint(ip)->setScatterer(nullptr); 
    trk->getPoint(ip)->deleteFitterInfo(rep);
    //TODO
    if (!trk->getPoint(ip)->hasRawMeasurements())
      trk->deletePoint(ip);
  }
}

void GblFitter::sortHits(Track* trk, const AbsTrackRep* rep) const { 
  // All measurement points in ref. track
  int npoints_meas = trk->getNumPointsWithMeasurement();  
  // Prepare state for extrapolation of track seed
  StateOnPlane reference(rep);
  rep->setTime(reference, trk->getTimeSeed());
  rep->setPosMom(reference, trk->getStateSeed());
  // Take the state to first plane
  SharedPlanePtr firstPlane(trk->getPointWithMeasurement(0)->getRawMeasurement(0)->constructPlane(reference));
  reference.extrapolateToPlane(firstPlane);
  //1st point is at arc-len=0
  double arcLenPos = 0;
  
  // Loop only between meas. points 
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas - 1; ipoint_meas++) {
    // current measurement point
    TrackPoint* point_meas = trk->getPointWithMeasurement(ipoint_meas);    
    // Current detector plane
    SharedPlanePtr plane = point_meas->getRawMeasurement(0)->constructPlane(reference);    
    // Get the next plane
    SharedPlanePtr nextPlane(trk->getPointWithMeasurement(ipoint_meas + 1)->getRawMeasurement(0)->constructPlane(reference));    
    
    point_meas->setSortingParameter(arcLenPos);
    arcLenPos += reference.extrapolateToPlane(nextPlane);
    
  } // end of loop over track points with measurement
  trk->getPointWithMeasurement(npoints_meas - 1)->setSortingParameter(arcLenPos);
  trk->sort();
}

std::vector<gbl::GblPoint> GblFitter::collectGblPoints(genfit::Track* trk, const genfit::AbsTrackRep* rep) {
  //TODO store collected points in in fit status? need streamer for GblPoint (or something like that)
  std::vector<gbl::GblPoint> thePoints;
  thePoints.clear();
  
  // Collect points from track and fitterInfo(rep)
  for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {   
    GblFitterInfo * gblfi = dynamic_cast<GblFitterInfo*>(trk->getPoint(ip)->getFitterInfo(rep));
    if (!gblfi)
      continue;
    thePoints.push_back(gblfi->constructGblPoint());      
  }  
  return thePoints;
}

void GblFitter::updateGblInfo(gbl::GblTrajectory& traj, genfit::Track* trk, const genfit::AbsTrackRep* rep) {
  //FIXME
  if (!traj.isValid())
    return;
  
  // Update points in track and fitterInfo(rep)
  int igblfi = -1;
  for (unsigned int ip = 0; ip < trk->getNumPoints(); ip++) {      
    GblFitterInfo * gblfi = dynamic_cast<GblFitterInfo*>(trk->getPoint(ip)->getFitterInfo(rep));
    if (!gblfi)
      continue;
    igblfi++;
    
    // The point will calculate its position on the track
    // (counting fitter infos) which hopefully
    gblfi->updateFitResults(traj);
    
    // This is agains logic. User can do this if he wants and it is recommended usually
    // so that following fit could reuse the updated seed
    //if (igblfi == 0) {
    //  trk->setStateSeed( gblfi->getFittedState(true).getPos(), gblfi->getFittedState(true).getMom() );
    //  trk->setCovSeed( gblfi->getFittedState(true).get6DCov() ); 
    //}    
  }
}

void GblFitter::getScattererFromMatList(double& length,
                             double& theta, double& s, double& ds,
                             const double p, const double mass, const double charge,
                             const std::vector<genfit::MatStep>& steps) const {
  theta = 0.; s = 0.; ds = 0.; length = 0;
  if (steps.empty()) return;
  
  // sum of step lengths
  double len = 0.;
  // normalization
  double sumxx = 0.;
  // first moment (non-normalized)
  double sumx2x2 = 0.;
  // (part of) second moment / variance (non-normalized)
  double sumx3x3 = 0.;
  
  double xmin = 0.;
  double xmax = 0.;
  
  for (unsigned int i = 0; i < steps.size(); i++) {
    const MatStep step = steps.at(i);
    // inverse of material radiation length ... (in 1/cm) ... "density of scattering"
    double rho = 1. / step.material_.radiationLength;
    len += fabs(step.stepSize_);
    xmin = xmax;
    xmax = xmin + fabs(step.stepSize_);
    // Compute integrals
    
    // integral of rho(x)
    sumxx   += rho * (xmax - xmin);
    // integral of x*rho(x)
    sumx2x2 += rho * (xmax * xmax - xmin * xmin) / 2.;
    // integral of x*x*rho(x)
    sumx3x3 += rho * (xmax * xmax * xmax - xmin * xmin * xmin) / 3.;
  }
  // This ensures PDG formula still gives positive results (but sumxx should be > 1e-4 for it to hold)
  if (sumxx < 1.0e-10) return;
  // Calculate theta from total sum of radiation length
  // instead of summimg squares of individual deflection angle variances
  // PDG formula:
  double beta = p / sqrt(p * p + mass * mass);
  theta = (0.0136 / p / beta) * fabs(charge) * sqrt(sumxx) * (1. + 0.038 * log(sumxx));
  //theta = (0.015 / p / beta) * fabs(charge) * sqrt(sumxx);
  
  // track length
  length = len;
  // Normalization factor
  double N = 1. / sumxx;
  // First moment
  s  = N * sumx2x2;
  // Square of second moment (variance)
  // integral of (x - s)*(x - s)*rho(x)
  double ds_2 = N * (sumx3x3 - 2. * sumx2x2 * s + sumxx * s * s);
  ds = sqrt(ds_2);
  
  #ifdef DEBUG
  ////std::cout << "Thick scatterer parameters (dtheta, <s>, ds): " << "(" << theta << ", " << s << ", " << ds << ")" << endl;  
  #endif
}

double GblFitter::constructGblInfo(Track* trk, const AbsTrackRep* rep)
{ 
  // All measurement points in ref. track
  int npoints_meas = trk->getNumPointsWithMeasurement();  
  // Dimesion of representation/state
  int dim = rep->getDim();  
  // Jacobian for point with measurement = how to propagate from previous point (scat/meas)
  TMatrixD jacPointToPoint(dim, dim);
  jacPointToPoint.UnitMatrix();
  
  // Prepare state for extrapolation of track seed
  // Take the state to first plane
  StateOnPlane reference(rep);
  rep->setTime(reference, trk->getTimeSeed());
  rep->setPosMom(reference, trk->getStateSeed());

  SharedPlanePtr firstPlane(trk->getPointWithMeasurement(0)->getRawMeasurement(0)->constructPlane(reference));  
  reference.extrapolateToPlane(firstPlane);
  
  double sumTrackLen = 0;
  // NOT used but useful
  TMatrixDSym noise; TVectorD deltaState;
  
  // Loop only between meas. points 
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    // current measurement point
    TrackPoint* point_meas = trk->getPointWithMeasurement(ipoint_meas);    
    // Current detector plane
    SharedPlanePtr plane = point_meas->getRawMeasurement(0)->constructPlane(reference);    
    // track direction at plane (in global coords)
    TVector3 trackDir = rep->getDir(reference);
    // track momentum direction vector at plane (in global coords)
    double trackMomMag = rep->getMomMag(reference);
    // charge of particle
    double particleCharge = rep->getCharge(reference);
    // mass of particle
    double particleMass = rep->getMass(reference);    
    // Parameters of a thick scatterer between measurements
    double trackLen = 0., scatTheta = 0., scatSMean = 0., scatDeltaS = 0.;
    // Parameters of two equivalent thin scatterers
    double theta1 = 0., theta2 = 0., s1 = 0., s2 = 0.;
    // jacobian from s1=0 to s2
    TMatrixD jacMeas2Scat(dim, dim);
    jacMeas2Scat.UnitMatrix();
    
    // Stop here if we are at last point (do not add scatterers to last point),
    // just the fitter info
    if (ipoint_meas >= npoints_meas - 1) {
            
      // Construction last measurement (no scatterer)
      // --------------------------------------------
      // Just add the fitter info of last plane
      GblFitterInfo* gblfimeas(new GblFitterInfo(point_meas, rep, reference));
      gblfimeas->setJacobian(jacPointToPoint);
      point_meas->setFitterInfo(gblfimeas);      
      // --------------------------------------------
      
      break;
    }
    // Extrapolate to next measurement to get material distribution
    // Use a temp copy of the StateOnPlane to propage between measurements
    StateOnPlane refCopy(reference);
    // Get the next plane
    SharedPlanePtr nextPlane(trk->getPointWithMeasurement(ipoint_meas + 1)->getRawMeasurement(0)->constructPlane(reference));    
    
    // Extrapolation for multiple scattering calculation
    // Extrap to point + 1, do NOT stop at boundary
    TVector3 segmentEntry = refCopy.getPos();
    rep->extrapolateToPlane(refCopy, nextPlane, false, false);
    TVector3 segmentExit = refCopy.getPos();
    
    getScattererFromMatList(trackLen,
                            scatTheta,
                            scatSMean,
                            scatDeltaS,
                            trackMomMag,
                            particleMass,
                            particleCharge,
                            rep->getSteps());
    // Now calculate positions and scattering variance for equivalent scatterers
    // (Solution from Claus Kleinwort (DESY))
    s1 = 0.; s2 = scatSMean + scatDeltaS * scatDeltaS / (scatSMean - s1);
    theta1 = sqrt(scatTheta * scatTheta * scatDeltaS * scatDeltaS / (scatDeltaS * scatDeltaS + (scatSMean - s1) * (scatSMean - s1)));
    theta2 = sqrt(scatTheta * scatTheta * (scatSMean - s1) * (scatSMean - s1) / (scatDeltaS * scatDeltaS + (scatSMean - s1) * (scatSMean - s1))); 
    
    // Call segment controller to set MS options:    
    if (m_segmentController)
      m_segmentController->controlTrackSegment(segmentEntry, segmentExit, this);    
    
    // Scattering options: OFF / THIN / THICK
    if (m_enableScatterers && !m_enableIntermediateScatterer) {
      theta1 = scatTheta;
      theta2 = 0;
    } else if (!m_enableScatterers) {
      theta1 = 0.;
      theta2 = 0.;
    }
    
    // Construction of measurement (with scatterer)
    // --------------------------------------------
    
    if (theta1 > scatEpsilon)  
      point_meas->setScatterer(new ThinScatterer(plane, Material(theta1, 0., 0., 0., 0.)));
    
    GblFitterInfo* gblfimeas(new GblFitterInfo(point_meas, rep, reference));
    gblfimeas->setJacobian(jacPointToPoint);
    point_meas->setFitterInfo(gblfimeas);
    // --------------------------------------------
    
    
    // If not last measurement, extrapolate and get jacobians for scattering points between this and next measurement
    if (theta2 > scatEpsilon) {
      // First scatterer has been placed at current measurement point (see above)      
      // theta2 > 0 ... we want second scatterer:
      // Extrapolate to s2 (we have s1 = 0)
      rep->extrapolateBy(reference, s2, false, true);  
      rep->getForwardJacobianAndNoise(jacMeas2Scat, noise, deltaState);
      
      // Construction of intermediate scatterer
      // --------------------------------------
      TrackPoint* scattp = new TrackPoint(trk);
      scattp->setSortingParameter(point_meas->getSortingParameter() + s2);
      scattp->setScatterer(new ThinScatterer(reference.getPlane(), Material(theta2, 0., 0., 0., 0.)));
      // Add point to track before next point
      int pointIndex = 0;
      //TODO Deduce this rather than looping over all points
      for (unsigned int itp = 0; itp < trk->getNumPoints(); itp++) {
        if (trk->getPoint(itp) == point_meas) {
          pointIndex = itp;
          break;
        }
      }
      trk->insertPoint(scattp, pointIndex + 1);      
      // Create and store fitter info
      GblFitterInfo * gblfiscat(new GblFitterInfo(scattp, rep, reference));
      gblfiscat->setJacobian(jacMeas2Scat);
      scattp->setFitterInfo(gblfiscat);
      // ---------------------------------------
            
      // Finish extrapolation to next measurement
      double nextStep = rep->extrapolateToPlane(reference, nextPlane, false, true);
      rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
      
      if (0. > nextStep) {
        cout << " ERROR: The extrapolation to measurement point " << (ipoint_meas + 2) << " stepped back by " << nextStep << "cm !!! Track will be cut before this point." << endl;
        // stop trajectory construction here
        break;
      }
      
    } else {
      // No scattering: extrapolate whole distance between measurements
      double nextStep = rep->extrapolateToPlane(reference, nextPlane, false, true);
      rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
      
      if (0. > nextStep) {
        cout << " ERROR: The extrapolation to measurement point " << (ipoint_meas + 2) << " stepped back by " << nextStep << "cm !!! Track will be cut before this point." << endl;
        // stop trajectory construction here
        break;
      }          
    }
    // Track length up to next point
    sumTrackLen += trackLen;    

  } // end of loop over track points with measurement
  return sumTrackLen;
}
