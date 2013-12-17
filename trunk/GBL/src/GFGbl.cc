//-*-mode: C++; c-basic-offset: 2; -*-
/* Copyright 2013
 *  Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *  This is an interface to General Broken Lines
 *
 *  Version: 3 (Tadeas)
 *  This version now supports both TrueHits and Clusters for VXD.
 *  It can be used for arbitrary material distribution between
 *  measurements. Moments of scattering distribution are computed
 *  and translated into two equivalent thin GBL scatterers placed
 *  at computed positions between measurement points.
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
// Avoid use of Kalman filter
#include "KalmanFitStatus.h"
#include "KalmanFittedStateOnPlane.h"

#include "MyDebugTools.h"
#include "Track.h"
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <list>
#include <FieldManager.h>
#include <HMatrixU.h>
#include <HMatrixV.h>


//#define DEBUG
#define OUTPUT
//#define KALMAN_CHECK


#ifdef DEBUG
ofstream debug("gbl.debug");
#endif

#ifdef OUTPUT
ofstream output("gbl.output");
#endif

// Millepede Binary File for output of GBL trajectories for alignment
gbl::MilleBinary* milleFile;

using namespace gbl;
using namespace std;
using namespace genfit;

GFGbl::GFGbl() :
  AbsFitter()
{
  milleFile = new MilleBinary("millefile.dat");
}

/**
 * @brief Evaluates moments of radiation length distribution from list of
 * material steps and computes parameters describing a corresponding thick scatterer.
 *
 * Based on input from Claus Kleinwort (DESY),
 * reformulated and adapted for continuous material distribution represented by
 * a sum of step functions. Returned thick scatterer can be represented by two GBL scattering points
 * at (s - ds) and (s + ds) with variance of theta equal to theta/sqrt(2) for both points.
 * Calculates variance of theta from total sum of radiation lengths
 * instead of summimg squares of individual deflection angle variances.
 *
 * @param theta returned: Variation of distribution of deflection angle
 * @param s returned: First moment of material scattering distribution
 * @param ds returned: Second moment (variance) of material scattering distribution
 * @param p Particle momentum magnitude (GeV/c)
 * @param mass Mass of particle (GeV/c)
 * @param steps Vector of material steps from (RKTrackRep) extrapolation
 * @return void
 */
void getScattererFromMatList(double& theta, double& s, double& ds, const double p, const double mass, const double charge, const std::vector<MatStep>& steps)
{
  theta = 0.; s = 0.; ds = 0.;
  if (steps.empty()) return;

  // normalization
  double sumxx = 0.;
  // first moment (non-normalized)
  double sumx2x2 = 0.;
  // (part of) second moment / variance (non-normalized)
  double sumx3x3 = 0.;

  double xmin = 0.;
  double xmax = 0.;

  for (MatStep step : steps) {
    // inverse of material radiation length (in 1/cm) ... "density of scattering"
    double radLen = 1. / step.materialProperties_.getRadLen();

    xmin = xmax;
    xmax = xmin + fabs(step.stepSize_);
    // Compute integrals
    sumxx   += radLen * (xmax - xmin);
    sumx2x2 += radLen * (xmax * xmax - xmin * xmin) / 2.;
    sumx3x3 += radLen * (xmax * xmax * xmax - xmin * xmin * xmin) / 3.;
  }
  // This ensures PDG formula still gives positive results (but sumxx should be > 1e-4 for it to hold)
  if (sumxx < 1.0e-10) return;
  // Calculate theta from total sum of radiation length
  // instead of summimg squares of individual deflection angle variances
  // PDG formula:
  double beta = p / sqrt(p * p + mass * mass);
  theta = (0.0136 / p / beta) * fabs(charge) * sqrt(sumxx) * (1. + 0.038 * log(sumxx));

  // Normalization factor
  double N = 1. / sumxx;
  // First moment
  s  = N * sumx2x2;
  // Second moment (variance)
  ds = N * (sumx3x3 - 2 * sumx2x2 * s + sumxx * s * s);

#ifdef DEBUG  
  cout << "Thick scatterer parameters:" << endl;
  cout << "Variance of theta: " << theta << endl;
  cout << "Mean s           : " << s << endl;
  cout << "Variance of s    : " << ds << endl;
#endif  
}

void GFGbl::processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits)
{
  // Chi2 of Reference Track
  double trkChi2 = 0.;

  // Dimesion of repository/state
  int dim = rep->getDim();
  // current measurement point
  TrackPoint* point_meas;
  // current raw measurement
  AbsMeasurement* raw_meas;

  int npoints_meas = trk->getNumPointsWithMeasurement();
  //  int npoints_all = trk->getNumPoints();
  
#ifdef DEBUG
  if (resortHits)
    cout << "WARNING: Hits resorting in GBL interface not supported." << endl;
  
  cout << "-------------------------------------------------------" << endl;
  cout << "               GBL processing genfit::Track            " << endl;
  cout << "-------------------------------------------------------" << endl;
  //  cout << " # Track Points       :  " << npoints_all  << endl;
  cout << " # Measurements Points:  " << npoints_meas << endl;
#endif
  
  std::vector<GblPoint> listOfPoints;
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

#ifdef KALMAN_CHECK
  //TODO: Remove this
  genfit::FitStatus* fs = trk->getFitStatus(rep);
  genfit::KalmanFitStatus* kfs = dynamic_cast<genfit::KalmanFitStatus*>(fs);
#endif
  
  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    point_meas = trk->getPointWithMeasurement(ipoint_meas);

    // Fitter info which contains Reference state and plane
    KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(point_meas->getFitterInfo(rep));
    if (!fi) {
      Exception e("KalmanFitterInfo (with reference state) for measurement does not exist", __LINE__, __FILE__);
      throw e;
    }
    // Current detector plane
    SharedPlanePtr plane = fi->getPlane();

    // StateOnPlane for extrapolation
    ReferenceStateOnPlane* reference = new ReferenceStateOnPlane(*fi->getReferenceState());
    // Representation state at plane
    TVectorD state = reference->getState();
    // track direction at plane (in global coords)
    TVector3 trackDir = rep->getDir(*reference);
    // track momentum vector at plane (in global coords)
    double trackMomMag = rep->getMomMag(*reference);
    // charge of particle
    double particleCharge = rep->getCharge(*reference);
    // mass of particle
    double particleMass = rep->getMass(*reference);

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
      AbsMeasurement* raw_measU;
      AbsMeasurement* raw_measV;

      int sensorId2 = 0;
      PlanarMeasurement* measPlanar2 = dynamic_cast<PlanarMeasurement*>(point_meas->getRawMeasurement(0));
      if (measPlanar2) sensorId2 = measPlanar->getPlaneId();

      // We only try to combine if at same sensor id (should be always, but who knows)
      // otherwise ignore this point
      if (sensorId != sensorId2) continue;

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
        continue;
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
    if (raw_meas->getRawHitCoords().GetNoElements() != 2) continue;

    // 2D hit coordinates
    TVectorD raw_coor = raw_meas->getRawHitCoords();
    // Covariance matrix of measurement
    TMatrixDSym raw_cov = raw_meas->getRawHitCov();
    // Projection matrix from repository state to measurement coords
    boost::scoped_ptr<const AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(rep));
    // Residual between measured position and reference track position
    TVectorD residual = raw_coor - HitHMatrix->Hv(state);

    trkChi2 += residual(0) * residual(0) / raw_cov(0, 0) + residual(1) * residual(1) / raw_cov(1, 1);

    double scatTheta = 0.;
    double scatSMean = 0.;
    double scatDeltaS = 0.;
    // Extrapolate to next measurement to get material distribution
    if (ipoint_meas < npoints_meas - 1) {
      // Extrap to point + 1, do NOT stop at boundary
      rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, false);
      getScattererFromMatList(scatTheta,
                              scatSMean,
                              scatDeltaS,
                              trackMomMag,
                              particleMass,
                              particleCharge,
                              rep->getSteps());
    }
    // Return back to state on current plane
    delete reference;
    reference = new ReferenceStateOnPlane(*fi->getReferenceState());

    TMatrixDSym noise;
    TVectorD deltaState;
    TMatrixD jacScat1ToScat2(dim, dim);
    TMatrixD jacScat2ToMeas(dim, dim);
    jacScat1ToScat2.UnitMatrix();
    jacScat2ToMeas.UnitMatrix();

    // If not last measurement, extrapolate and get jacobians for scattering points between this and next measurement
    if (ipoint_meas < npoints_meas - 1) {
      // Only add scatteres if theta>0
      if (scatTheta > 0.) {
        // Extrapolate to (s - ds)
        rep->extrapolateBy(*reference, scatSMean - scatDeltaS, false);
        rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
        // Extrapolate to (s + ds)
        rep->extrapolateBy(*reference, 2.0 * scatDeltaS, false);
        rep->getForwardJacobianAndNoise(jacScat1ToScat2, noise, deltaState);
        // Finish extrapolation to next measurement
        rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, true);
        rep->getForwardJacobianAndNoise(jacScat2ToMeas, noise, deltaState);
      } else {
        // No scattering: extrapolate whole distance between measurements
        rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas + 1)->getFitterInfo(rep)->getPlane(), false, true);
        rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
      }
    }
    // Measurement point
    GblPoint measPoint(jacPointToPoint);
    // Projection from local (state) coordinates to measurement coordinates (inverted)
    TMatrixD proL2m(2, 2);
    proL2m(0, 0) = HitHMatrix->getMatrix()(0, 3);
    proL2m(0, 1) = HitHMatrix->getMatrix()(0, 4);
    proL2m(1, 0) = HitHMatrix->getMatrix()(1, 3);
    proL2m(1, 1) = HitHMatrix->getMatrix()(1, 4);
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
    n_gbl_points++;
    n_gbl_meas_points++;

    // Now store scatterers if not last measurement
    if (scatTheta > 0. && ipoint_meas < npoints_meas - 1) {
      // TrackRep state is perpendicular to track direction if using extrapolateBy (I asked Johannes Rauch),
      // therefore scattering covariance is diagonal and and both elements are equal
      TMatrixDSym scatCov(2);
      scatCov.Zero();
      scatCov(0, 0) = scatTheta * scatTheta / 2.0;
      scatCov(1, 1) = scatTheta * scatTheta / 2.0;
      // Now invert scattering covariance to get scattering precision
      scatCov.Invert();
      // We assume no initial kinks
      TVectorD scatResidual(2);
      scatResidual.Zero();

      GblPoint scatPoint1(jacScat1ToScat2);
      scatPoint1.addScatterer(scatResidual, scatCov);
      listOfPoints.push_back(scatPoint1);
      n_gbl_points++;

      GblPoint scatPoint2(jacScat2ToMeas);
      scatPoint2.addScatterer(scatResidual, scatCov);
      listOfPoints.push_back(scatPoint2);
      n_gbl_points++;
    }
    // Free memory on the heap
    delete reference;
  }
  // We should have at least two measurement points to fit anything
  if (n_gbl_meas_points > 1) {
    GblTrajectory* traj = 0;
    try {
      //TODO: Use clever way to determine zero B-field
      double Bfield = genfit::FieldManager::getInstance()->getFieldVal(TVector3(0., 0., 0.)).Mag();
      bool fitQoverP = Bfield > 0.;
      
      traj = new GblTrajectory(listOfPoints, seedLabel, clSeed, fitQoverP);
      traj->fit(Chi2, Ndf, lostWeight);
    } catch (...) {
      // Gbl failed critically (usually GblPoint::getDerivatives ... singular matrix inversion)
      delete traj;
      return;
    }
    // GBL fit succeded if Ndf > 0
    //TODO: Here hould be some track quality check
    if (Ndf > 1) {
      traj->milleOut(*milleFile);
      delete traj;
      
      //TODO Second GBL external iteration ... for electrons with Brehmstrahlung
      
#ifdef DEBUG
      cout << "____________GBL Fit Results________________" << endl;
      cout << "Initial Chi2:" << trkChi2 << endl;
      cout << "Chi2:        " << Chi2 << endl;
      cout << "NDF:         " << Ndf << endl;      
#endif

#ifdef OUTPUT
  #ifdef KALMAN_CHECK
      //TODO: Avoid use of Kalman filter
      output << kfs->getBackwardNdf() << " " << kfs->getBackwardChi2() << " ";
  #else
      output << 1 << " " << 0. << " " << endl;          
  #endif      

      output << Ndf << " " << Chi2 << " " << endl;    
#endif
    }
    
  }
}
