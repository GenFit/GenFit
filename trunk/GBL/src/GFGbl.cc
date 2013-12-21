//-*-mode: C++; c-basic-offset: 2; -*-
/* Copyright 2013
 *  Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *  This is an interface to General Broken Lines
 *
 *  Version: 4 (Tadeas)
 *  Fixed calculation of equvivalent scatterers (solution by C. Kleinwort)
 *  Now a scatterer is inserted at each measurement (except last) and between each two measurements.
 *  TrueHits/Clusters. Ghost (1D) hits ignored. With or without magnetic field.
 *  Version: 3 (Tadeas)
 *  This version now supports both TrueHits and Clusters for VXD.
 *  It can be used for arbitrary material distribution between
 *  measurements. Moments of scattering distribution are computed
 *  and translated into two equivalent thin GBL scatterers placed
 *  at computed positions between measurement points.
 *  Version: 2 ... never published (Tadeas)
 *  Scatterer at each boundary (tooo many scatterers). TrueHits/Clusters. Without global der.&MP2 output
 *  Version: 1 (Sergey & Tadeas)
 *  Scatterers at measurement planes. TrueHits
 *  Version 0: (Sergey)
 *  Without scatterers. Genfit 1.
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

#include "GFGbl.h"
#include "GblTrajectory.h"
#include "GblPoint.h"

#include "AbsMeasurement.h"
#include "PlanarMeasurement.h"
#include "KalmanFitterInfo.h"

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


#define DEBUG
#define OUTPUT


#ifdef DEBUG
ofstream debug("gbl.debug");
TFile* diag;
TH1F* resHistosU[6];
TH1F* resHistosV[6];
TH1F* mhistosU[6];
TH1F* mhistosV[6];
TH1F* ghistosU[6];
TH1F* ghistosV[6];
TH1F* chi2OndfHisto;
TH1F* pValueHisto;


#endif

#ifdef OUTPUT
void writeHistoDataForLabel(double label, TVectorD res, TVectorD measErr, TVectorD resErr, TVectorD downWeights) {
  if (label > 6) return;
  if (label < 1) return;
  int i = int(label);
  resHistosU[i-1]->Fill(res[0]);
  resHistosV[i-1]->Fill(res[1]);
  mhistosU[i-1]->Fill(res[0]/measErr[0]);
  mhistosV[i-1]->Fill(res[1]/measErr[1]);
  ghistosU[i-1]->Fill(res[0]/resErr[0]);
  ghistosV[i-1]->Fill(res[1]/resErr[1]);        
}
#endif

// Millepede Binary File for output of GBL trajectories for alignment
gbl::MilleBinary* milleFile;
// Minimum scattering sigma (will be inverted...)
const double scatEpsilon = 1.e-8;


using namespace gbl;
using namespace std;
using namespace genfit;

GFGbl::GFGbl() :
AbsFitter()
{
  milleFile = new MilleBinary("millefile.dat");
  
#ifdef OUTPUT
  diag = new TFile("gbl.root","RECREATE"); 
  char name[20];
  
  for (int i = 0; i < 6; i++){
    sprintf(name, "res_u_%i", i+1);
    resHistosU[i] = new TH1F(name, "Residual (U)", 1000, -0.1, 0.1);
    sprintf(name, "res_v_%i", i+1);
    resHistosV[i] = new TH1F(name, "Residual (V)", 1000, -0.1, 0.1);
    sprintf(name, "meas_pull_u_%i", i+1);
    mhistosU[i] = new TH1F(name, "Res/Meas.Err. (U)", 1000, -10., 10.);
    sprintf(name, "meas_pull_v_%i", i+1);
    mhistosV[i] = new TH1F(name, "Res/Meas.Err. (V)", 1000, -10., 10.);
    sprintf(name, "pull_u_%i", i+1);
    ghistosU[i] = new TH1F(name, "Res/Res.Err. (U)", 1000, -10., 10.);
    sprintf(name, "pull_v_%i", i+1);
    ghistosV[i] = new TH1F(name, "Res/Res.Err. (V)", 1000, -10., 10.);
  }
  chi2OndfHisto = new TH1F("chi2_ndf", "Track Chi2/NDF", 100, 0., 10.);
  pValueHisto = new TH1F("p_value", "Track P-value", 100, 0., 1.);
#endif  
}

void GFGbl::endRun()
{
#ifdef OUTPUT
  diag->cd();
  diag->Write();
  diag->Close();
#endif
  // This is needed to close the file before alignment starts
  delete milleFile;
}

/**
 * @brief Evaluates moments of radiation length distribution from list of
 * material steps and computes parameters describing a corresponding thick scatterer.
 *
 * Based on input from Claus Kleinwort (DESY),
 * adapted for continuous material distribution represented by
 * a sum of step functions. Returned thick scatterer can be represented by two GBL scattering points
 * at (s - ds) and (s + ds) with variance of theta equal to theta/sqrt(2) for both points.
 * Calculates variance of theta from total sum of radiation lengths
 * instead of summimg squares of individual deflection angle variances.
 *
 * @param length returned: Length of the track
 * @param theta returned: Variation of distribution of deflection angle
 * @param s returned: First moment of material scattering distribution
 * @param ds returned: Second moment (variance) of material scattering distribution
 * @param p Particle momentum magnitude (GeV/c)
 * @param mass Mass of particle (GeV/c/c)
 * @param steps Vector of material steps from (RKTrackRep) extrapolation
 * @return void
 */
void getScattererFromMatList(double& length, double& theta, double& s, double& ds, const double p, const double mass, const double charge, const std::vector<MatStep>& steps)
{
  theta = 0.; s = 0.; ds = 0.;
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
    double rho = 1. / step.materialProperties_.getRadLen();
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
  //cout << "Thick scatterer parameters:" << endl;
  //cout << "Variance of theta: " << theta << endl;
  //cout << "Mean s           : " << s << endl;
  //cout << "Variance of s    : " << ds << endl;
  
  #endif
}

void GFGbl::processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits)
{
  // Flag used to mark error in raw measurement combination
  // measurement won't be considered, but scattering yes
  bool skipMeasurement = false;
  // Chi2 of Reference Track
  double trkChi2 = 0.;
  // This flag enables/disables fitting of q/p parameter in GBL
  // It is switched off automatically if no B-field at (0,0,0) is detected.
  bool fitQoverP = true;
  //TODO: Use clever way to determine zero B-field
  double Bfield = genfit::FieldManager::getInstance()->getFieldVal(TVector3(0., 0., 0.)).Mag();
  if (!(Bfield > 0.))
    fitQoverP = false;
  
  // Dimesion of repository/state
  int dim = rep->getDim();
  // current measurement point
  TrackPoint* point_meas;
  // current raw measurement
  AbsMeasurement* raw_meas;
  
  // We assume no initial kinks, this will be reused several times
  TVectorD scatResidual(2);
  scatResidual.Zero();
  
  // All measurement points in ref. track
  int npoints_meas = trk->getNumPointsWithMeasurement();
  
  #ifdef DEBUG
  int npoints_all = trk->getNumPoints();
  
  if (resortHits)
    cout << "WARNING: Hits resorting in GBL interface not supported." << endl;
  
  cout << "-------------------------------------------------------" << endl;
  cout << "               GBL processing genfit::Track            " << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << " # Ref. Track Points  :  " << npoints_all  << endl;
  cout << " # Meas. Points       :  " << npoints_meas << endl;
  
  #endif
  // List of prepared GBL points for GBL trajectory construction
  std::vector<GblPoint> listOfPoints;
  
  std::vector<double> listOfLayers;
  //TODO: Add internal/external seed (from CDC) option in the future
  // index of point with seed information (0 for none)
  unsigned int seedLabel = 0;
  // Seed covariance
  // TMatrixDSym clCov(dim);
  // Seed state
  TMatrixDSym clSeed(dim);
  
  // propagation Jacobian to next point from current measurement point
  TMatrixD jacPointToPoint(dim, dim);
  jacPointToPoint.UnitMatrix();
  
  int n_gbl_points = 0;
  int n_gbl_meas_points = 0;
  int Ndf = 0;
  double Chi2 = 0.;
  double lostWeight = 0.;
  
  // Momentum of track at current plane
  double trackMomMag = 0.;
  // Charge of particle at current plane :-)
  double particleCharge = 1.;
  
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    point_meas = trk->getPointWithMeasurement(ipoint_meas);
    
    if (!point_meas->hasFitterInfo(rep)) {
      cout << " ERROR: Measurement point does not have a fitter info. Track will be skipped." << endl;
      return;
    }
    // Fitter info which contains Reference state and plane
    KalmanFitterInfo* fi = dynamic_cast<KalmanFitterInfo*>(point_meas->getFitterInfo(rep));
    if (!fi) {
      cout << " ERROR: KalmanFitterInfo (with reference state) for measurement does not exist. Track will be skipped." << endl;
      return;
    }
    // Current detector plane
    SharedPlanePtr plane = fi->getPlane();    
    if (!fi->hasReferenceState()) {
      cout << " ERROR: Fitter info does not contain reference state. Track will be skipped." << endl;
      return;
    }
    // Reference StateOnPlane for extrapolation
    ReferenceStateOnPlane* reference = new ReferenceStateOnPlane(*fi->getReferenceState());
    // Representation state at plane
    TVectorD state = reference->getState();
    // track direction at plane (in global coords)
    TVector3 trackDir = rep->getDir(*reference);
    // track momentum vector at plane (in global coords)
    trackMomMag = rep->getMomMag(*reference);
    // charge of particle
    particleCharge = rep->getCharge(*reference);
    // mass of particle
    double particleMass = rep->getMass(*reference);
    
    // Parameters of a thick scatterer between measurements
    double trackLen = 0.;
    double scatTheta = 0.;
    double scatSMean = 0.;
    double scatDeltaS = 0.;
    // Parameters of two equivalent thin scatterers
    double theta1 = 0.;
    double theta2 = 0.;
    double s1 = 0.;
    double s2 = 0.;
    
    TMatrixDSym noise;
    TVectorD deltaState;
    // jacobian from s2 to M2
    TMatrixD jacScat2ToMeas(dim, dim);
    jacScat2ToMeas.UnitMatrix();
    
    //cout << " Starting extrapolation..." << endl;
    try {

    // Extrapolate to next measurement to get material distribution
    if (ipoint_meas < npoints_meas - 1) {
      // Check if fitter info is in place
      if (!trk->getPoint(ipoint_meas + 1)->hasFitterInfo(rep)) {
        cout << " ERROR: Measurement point does not have a fitter info. Track will be skipped." << endl;
        return;
      }
      // Fitter of next point info which is only used now to get the plane
      KalmanFitterInfo* fi_i_plus_1 = dynamic_cast<KalmanFitterInfo*>(trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep));
      if (!fi_i_plus_1) {
        cout << " ERROR: KalmanFitterInfo (with reference state) for measurement does not exist. Track will be skipped." << endl;
        return;
      }
      
      // Extrap to point + 1, do NOT stop at boundary
      rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, false);
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
      s1 = 0.;
      s2 = scatSMean + scatDeltaS * scatDeltaS / (scatSMean - s1);
      theta1 = sqrt(scatTheta * scatTheta * scatDeltaS * scatDeltaS / (scatDeltaS * scatDeltaS + (scatSMean - s1) * (scatSMean - s1)));
      theta2 = sqrt(scatTheta * scatTheta * (scatSMean - s1) * (scatSMean - s1) / (scatDeltaS * scatDeltaS + (scatSMean - s1) * (scatSMean - s1)));
      
      if (s2 >= trackLen - 1.e-16 || s2 <= 1.e-16) {
        cout << " ERROR: Failed GBL point position calculation. Points too close. GBLTrajectory construction might fail. Let's try it..." << endl;
      }
      
    }
    // Return back to state on current plane
    delete reference;
    reference = new ReferenceStateOnPlane(*fi->getReferenceState());
    
    // If not last measurement, extrapolate and get jacobians for scattering points between this and next measurement
    if (ipoint_meas < npoints_meas - 1) {
      if (theta2 > scatEpsilon) {
        // First scatterer will be placed at current measurement point (see bellow)
        
        // theta2 > 0 ... we want second scatterer:
        // Extrapolate to s2 (remember s1 = 0)
        rep->extrapolateBy(*reference, s2, false);
        rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
        // Finish extrapolation to next measurement
        //TODO: check if fitter info at next point exists...
        rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, true);
        rep->getForwardJacobianAndNoise(jacScat2ToMeas, noise, deltaState);
      } else {
        // No scattering: extrapolate whole distance between measurements
        rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, true);
        rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
      }
    }
    } catch(...) {
      cout << " ERROR: Extrapolation failed. Track will be skipped." << endl;
      return;
    }
    //cout << " Extrapolation finished." << endl;
    
    
    // Now get measurement. First have a look if we need to combine SVD clusters...
    
    // Try to get VxdId of current plane
    int sensorId = 0;
    PlanarMeasurement* measPlanar = dynamic_cast<PlanarMeasurement*>(point_meas->getRawMeasurement(0));
    if (measPlanar) sensorId = measPlanar->getPlaneId();
    
    //WARNING: Now we only support 2D measurements. If 2 raw measurements are stored at the track
    // point, these are checked if they correspond to "u" and "v" measurement (for SVD clusters) and these
    // measurements are combined. SLANTED SENSORS NOT YET SUPPORTED!!!
    if (point_meas->getRawMeasurement(0)->getDim() != 2
      && trk->getPointWithMeasurement(ipoint_meas)->getNumRawMeasurements() == 2
      && point_meas->getRawMeasurement(0)->getDim() == 1
      && point_meas->getRawMeasurement(1)->getDim() == 1) {
    AbsMeasurement* raw_measU = 0;
    AbsMeasurement* raw_measV = 0;
    
    // cout << " Two 1D Measurements encountered. " << endl;
    
    int sensorId2 = -1;
    PlanarMeasurement* measPlanar2 = dynamic_cast<PlanarMeasurement*>(point_meas->getRawMeasurement(0));
    if (measPlanar2) sensorId2 = measPlanar->getPlaneId();
    
    // We only try to combine if at same sensor id (should be always, but who knows)
    // otherwise ignore this point
    if (sensorId != sensorId2) {
      skipMeasurement = true;
      cout << " ERROR: Incompatible sensorIDs at measurement point " << ipoint_meas << ". Measurement will be skipped." << endl;
    }
    
    // We have to combine two SVD 1D Clusters at the same plane into one 2D recohit
    AbsMeasurement* raw_meas1 = point_meas->getRawMeasurement(0);
    AbsMeasurement* raw_meas2 = point_meas->getRawMeasurement(1);
    // Decide which cluster is u and which v based on H-matrix
    if (raw_meas1->constructHMatrix(rep)->isEqual(genfit::HMatrixU())
      && raw_meas2->constructHMatrix(rep)->isEqual(genfit::HMatrixV())) {
      // right order U, V
      raw_measU = raw_meas1;
    raw_measV = raw_meas2;
      } else if (raw_meas2->constructHMatrix(rep)->isEqual(genfit::HMatrixU())
        && raw_meas1->constructHMatrix(rep)->isEqual(genfit::HMatrixV())) {
        // inversed order V, U
        raw_measU = raw_meas2;
      raw_measV = raw_meas1;
        } else {
          // Incompatible measurements ... skip this point
          skipMeasurement = true;
          cout << " ERROR: Incompatible 1D measurements at meas. point " << ipoint_meas << ". Measurement will be skipped." << endl;
          
        }
        // Combine raw measurements
        TVectorD _raw_coor(2);
        _raw_coor(0) = raw_measU->getRawHitCoords()(0);
        _raw_coor(1) = raw_measV->getRawHitCoords()(0);
        // Combine covariance matrix
        TMatrixDSym _raw_cov(2);
        _raw_cov.Zero();
        _raw_cov(0, 0) = raw_measU->getRawHitCov()(0, 0);
        _raw_cov(1, 1) = raw_measV->getRawHitCov()(0, 0);
        // Create new combined measurement
        raw_meas = new PlanarMeasurement(_raw_coor, _raw_cov, raw_measU->getDetId(), raw_measU->getHitId(), point_meas);
      } else {
        // Default behavior
        raw_meas = point_meas->getRawMeasurement(0);
      }
      //TODO: We only support 2D measurements in GBL (ot two 1D combined above)
      if (raw_meas->getRawHitCoords().GetNoElements() != 2) {
        skipMeasurement = true;
        #ifdef DEBUG
        cout << " WARNING: Measurement " << (ipoint_meas + 1) << " is not 2D. Will be skipped. " << endl;
        #endif
      }
      
      // Now we have all necessary information, so lets insert current measurement point
      // if we don't want to skip it (e.g. ghost SVD hit ... just 1D information)
      if (!skipMeasurement) {
        // 2D hit coordinates
        TVectorD raw_coor = raw_meas->getRawHitCoords();
        // Covariance matrix of measurement
        TMatrixDSym raw_cov = raw_meas->getRawHitCov();
        // Projection matrix from repository state to measurement coords
        boost::scoped_ptr<const AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(rep));
        // Residual between measured position and reference track position
        TVectorD residual = raw_coor - HitHMatrix->Hv(state);
        
        trkChi2 += residual(0) * residual(0) / raw_cov(0, 0) + residual(1) * residual(1) / raw_cov(1, 1);
        
        // Measurement point
        GblPoint measPoint(jacPointToPoint);
        // Projection from local (state) coordinates to measurement coordinates (inverted)
        // 2x2 matrix ... last block of H matrix (2 rows x 5 columns)
        TMatrixD proL2m = HitHMatrix->getMatrix().GetSub(0, 1, 3, 4);
        proL2m.Invert();
        measPoint.addMeasurement(proL2m, residual, raw_cov.Invert());
        
        //Add global derivatives to the point
        
        // sensor label = sensorID * 10, then last digit is label for global derivative for the sensor
        int label = sensorId * 10;
        // values for global derivatives
        TMatrixD derGlobal(2, 6);
        // labels for global derivatives
        std::vector<int> labGlobal;
        
        // track direction in global coords
        TVector3 tDir = trackDir;
        // sensor u direction in global coords
        TVector3 uDir = plane->getU();
        // sensor v direction in global coords
        TVector3 vDir = plane->getV();
        // sensor normal direction in global coords
        TVector3 nDir = plane->getNormal();
        // track direction in local sensor system
        TVector3 tLoc = TVector3(uDir.Dot(tDir), vDir.Dot(tDir), nDir.Dot(tDir));
        
        // track u-slope in local sensor system
        double uSlope = tLoc[0] / tLoc[2];
        // track v-slope in local sensor system
        double vSlope = tLoc[1] / tLoc[2];
        
        // Measured track u-position in local sensor system
        double uPos = raw_coor[0];
        // Measured track v-position in local sensor system
        double vPos = raw_coor[1];
        
        //Global derivatives for alignment in sensor local coordinates
        derGlobal(0, 0) = 1.0;
        derGlobal(0, 1) = 0.0;
        derGlobal(0, 2) = - uSlope;
        derGlobal(0, 3) = vPos * uSlope;
        derGlobal(0, 4) = -uPos * uSlope;
        derGlobal(0, 5) = vPos;
        
        derGlobal(1, 0) = 0.0;
        derGlobal(1, 1) = 1.0;
        derGlobal(1, 2) = - vSlope;
        derGlobal(1, 3) = vPos * vSlope;
        derGlobal(1, 4) = -uPos * vSlope;
        derGlobal(1, 5) = -uPos;
        
        labGlobal.push_back(label + 1); // u
        labGlobal.push_back(label + 2); // v
        labGlobal.push_back(label + 3); // w
        labGlobal.push_back(label + 4); // alpha
        labGlobal.push_back(label + 5); // beta
        labGlobal.push_back(label + 6); // gamma
        
        measPoint.addGlobals(labGlobal, derGlobal);
        
        listOfPoints.push_back(measPoint);
        listOfLayers.push_back((unsigned int) sensorId >> 13);
        n_gbl_points++;
        n_gbl_meas_points++;
      } else {
        // Incompatible measurement, store point without measurement
        GblPoint dummyPoint(jacPointToPoint);
        listOfPoints.push_back(dummyPoint);
        listOfLayers.push_back((unsigned int) sensorId >> 13);
        n_gbl_points++;
        skipMeasurement = false;
        #ifdef DEBUG
        cout << " Dummy point inserted. " << endl;
        #endif
      }
      
      // Now store scatterers if not last measurement and if we decided
      // there should be scatteres, otherwise the jacobian in measurement
      // stored above is already correct
      if (ipoint_meas < npoints_meas - 1) {
        
        if (theta1 > scatEpsilon) {
          // We have to insert first scatterer at measurement point
          // Therefore (because state is perpendicular to plane, NOT track)
          // we have non-diaonal matrix of multiple scattering covariance
          // We have to project scattering into plane coordinates
          double c1 = trackDir.Dot(plane->getU());
          double c2 = trackDir.Dot(plane->getV());
          TMatrixDSym scatCov(2);
          scatCov(0, 0) = 1. - c2 * c2;
          scatCov(1, 1) = 1. - c1 * c1;
          scatCov(0, 1) = c1 * c2;
          scatCov(1, 0) = c1 * c2;
          scatCov *= theta1 * theta1 / (1. - c1 * c1 - c2 * c2) / (1. - c1 * c1 - c2 * c2) ;
          
          // last point is the just inserted measurement (or dummy point)
          GblPoint& lastPoint = listOfPoints.at(ipoint_meas);
          lastPoint.addScatterer(scatResidual, scatCov.Invert());
          
        }
        
        if (theta2 > scatEpsilon) {
          // second scatterer is somewhere between measurements
          // TrackRep state is perpendicular to track direction if using extrapolateBy (I asked Johannes Rauch),
          // therefore scattering covariance is diagonal (and both elements are equal)
          TMatrixDSym scatCov(2);
          scatCov.Zero();
          scatCov(0, 0) = theta2 * theta2;
          scatCov(1, 1) = theta2 * theta2;
          
          GblPoint scatPoint(jacScat2ToMeas);
          scatPoint.addScatterer(scatResidual, scatCov.Invert());
          listOfPoints.push_back(scatPoint);
          listOfLayers.push_back(((unsigned int) sensorId >> 13)+0.5);
          n_gbl_points++;
        }
        
        
      }
      // Free memory on the heap
      delete reference;
  }
  // We should have at least two measurement points to fit anything
  if (n_gbl_meas_points > 1) {
    bool writeOut = true;
    int fitRes = 0;
    GblTrajectory* traj = 0;
    try {
      // Construct the GBL trajectory, seed not used
      traj = new GblTrajectory(listOfPoints, seedLabel, clSeed, fitQoverP);
      // Fit the trajectory
      fitRes = traj->fit(Chi2, Ndf, lostWeight);
      if (fitRes != 0)
        writeOut = false;
    } catch (...) {
      // Gbl failed critically (usually GblPoint::getDerivatives ... singular matrix inversion)
      writeOut = false;
    }
    
#ifdef DEBUG
    cout << " Ref. Track Chi2      :  " << trkChi2 << endl;
    cout << " Ref. end momentum    :  " << trackMomMag << " GeV/c ";
    if (abs(trk->getCardinalRep()->getPDG()) == 11) {
      if (particleCharge < 0.)
        cout << "(electron)";
      else
        cout << "(positron)";
    }
    cout << endl;
    
    cout << "------------------ GBL Fit Results --------------------" << endl;
    cout << " Fit q/p parameter    :  " << ((fitQoverP) ? ("True") : ("False")) << endl;
    cout << " Valid trajectory     :  " << ((traj->isValid()) ? ("True") : ("False")) << endl;
    cout << " Fit result           :  " << fitRes << "    (0 for success)" << endl;
    cout << " # GBL meas. points   :  " << n_gbl_meas_points << endl;
    cout << " # GBL all points     :  " << n_gbl_points << endl;
    cout << " GBL track NDF        :  " << Ndf << "    (-1 for failure)" << endl;
    cout << " GBL track Chi2       :  " << Chi2 << endl;
    cout << " GBL track P-value    :  " << TMath::Prob(Chi2, Ndf) << endl;
    cout << "-------------------------------------------------------" << endl;
    //traj->printTrajectory(100);
    //traj->printData();
    //traj->printPoints(100);      
#endif
    
    // GBL fit succeded if Ndf > 0
    //TODO: Here hould be some track quality check
    if (Ndf > 1 && traj->isValid() && writeOut) {
      // Write trajectory data to Millepede II binary file
      traj->milleOut(*milleFile);
      
#ifdef DEBUG
      cout << " GBL Track written to Millepede II binary file." << endl;
      cout << "-------------------------------------------------------" << endl;      
#endif
      
#ifdef OUTPUT
      // Fill histograms
      chi2OndfHisto->Fill(Chi2/Ndf);
      pValueHisto->Fill(TMath::Prob(Chi2, Ndf));
      // Loop over all GBL points
      for (unsigned int p = 0; p < listOfPoints.size(); p++){
        unsigned int label = p+1;
        unsigned int numRes;
        TVectorD residuals(2);
        TVectorD measErrors(2);
        TVectorD resErrors(2);
        TVectorD downWeights(2);
        //TODO: now we only provide info about measurements, not kinks
        if (!listOfPoints.at(p).hasMeasurement())
          continue;
        
        traj->getMeasResults(label, numRes, residuals, measErrors, resErrors, downWeights);
        // Write layer-wise data 
        //TODO: we only support 6 layers now
        writeHistoDataForLabel(listOfLayers.at(p), residuals, measErrors, resErrors, downWeights);        
      } // end for points
      
#endif
    } // end if "all succeded"
    
    // Free memory
    delete traj;    
  }
}


