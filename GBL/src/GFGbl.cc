//-*-mode: C++; c-basic-offset: 2; -*-
/* Copyright 2013
   Authors: Sergey Yashchenko and Tadeas Bilka
 
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

#include "GFGbl.h"
#include "GblTrajectory.h"
#include "GblPoint.h"

#include "AbsMeasurement.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"
#include "KalmanFittedStateOnPlane.h"
#include "MyDebugTools.h"
#include "Track.h"

#include <TGeoManager.h>
#include <string>
#include <list>
#include <FieldManager.h>

#define DEBUG
#define OUTPUT

#ifdef DEBUG
ofstream ofs_debug("gbl.debug");
#endif

#ifdef OUTPUT
ofstream ofs_output("gbl.output");
#endif

// Millepede Binary File for output of GBL trajectories for alignment
gbl::MilleBinary *milleFile;

using namespace gbl;
using namespace std;
using namespace genfit;

GFGbl::GFGbl() :
  AbsFitter()
{
  milleFile = new MilleBinary("vxd.dat");
}

/**
 * @brief Simple Jacobian in the limit of q/p -> 0 for magnetic field B = (0, 0, Bz)
 *
 * @param ds Length of path of the particle passed during extrapolation
 * @param Bz Magnetic field component in z-coord
 * @param cosLambda Cosinus of "lambda" angle ... got from getCosLambda(particle momentum at plane where extrap. starts)
 * @return TMatrixD
 */

TMatrixD getSimpleJac(double ds, double Bz, double cosLambda) {
  TMatrixD j(5, 5);
  j.UnitMatrix();

  j(1, 0) = - Bz * ds * cosLambda;
  j(3, 0) = - Bz * 0.5 * ds * ds;
  j(3, 1) = ds;
  j(4, 2) = ds;
  return j;
}

/**
 * @brief Extracts basf2's VXD GeoCache sensorId from geometry node name
 *
 * @param path Path of the node of active volume of the sensor
 * @return int
 */
int getSensorIdFromPath(const char *path) {
  // Inspired by http://stackoverflow.com/questions/236129/splitting-a-string-in-c
  // Answer #562
  cout << path << endl;
  std::string l(path);
  std::replace(l.begin(), l.end(), '_', ' ');
  std::istringstream stm(l);
  std::vector<std::string> tokens;
  for (;;) {
    std::string word;
    if (!(stm >> word))
      break;
    tokens.push_back(word);
  }
  int res = atoi(tokens[tokens.size()-1].c_str());
  if (res < 1000 && tokens.size()>=2)
    res = atoi(tokens[tokens.size()-2].c_str());
  return res;
}

/**
 * @brief Returns squared variance of the scattering angle
 *
 * @param momMag Particle's momentum magnitude [GeV]
 * @param XX0 Fraction of radiation length (x/X0) passed in sensor material
 * @return double
 */
double getScattSigma2(double momMag, double XX0) {

  double momFactor = (0.0136 / momMag);
  double x0Factor = sqrt(XX0) * (1. + 0.038 * log(XX0));

  double scatSigma = momFactor * x0Factor;

  return scatSigma * scatSigma;
}

/**
 * @brief Returns sinus of angle "phi" of track projection into x-y (R-Phi) plane
 *
 * @param trackMom ...
 * @return double
 */
double getSinPhi(TVector3 trackMom) {
  return sin(atan(trackMom[1]/trackMom[0]));
}

/**
 * @brief Returns cosinus of the angle "lambda" of track projection into z-plane
 *
 * @param trackMom Vector of track momentum
 * @return double
 */
double getCosLambda(TVector3 trackMom) {
  double cosLambda = 1.;
  double pt = sqrt(trackMom[0] * trackMom[0] + trackMom[1] * trackMom[1]);
  double dzds = trackMom[2] / pt;
  cosLambda = 1. / sqrt(1. + dzds * dzds);
  return cosLambda;
}

/**
 * @brief Returns squared cosinus of the angle "lambda" of track projection into z-plane
 *
 * @param trackMom Vector of track momentum
 * @return double
 */
double getCosLambda2(TVector3 trackMom) {
  double cosLambda = getCosLambda(trackMom);
  return cosLambda*cosLambda;
}

/**
 * @brief Jacobian for transformation from local to curvilinear (co-moving) frame. B field is retrived at sensor plane
 * See paper from Strandlie & Wittek
 * Nuclear Instruments and Methods in Physics Research
 * Derivation of Jacobians for the propagation of covariance matrices of track parameters in homogeneous magnetic fields
 * A 566 (2006) 687-698
 * page 697, formula A.24
 *
 * @param plane SharedPlanePtr for extraction of normal & measurement directions
 * @param mom TVector3 vector of particle momentum at sensor plane
 * @param charge +1 or -1 ... charge of particle
 * @return TMatrixD
 */
TMatrixD jacLocal2Curvilinear(SharedPlanePtr plane, TVector3 mom, int charge, TVector3 pos) {
  TVector3 B = 0.3 / 100. * genfit::FieldManager::getInstance()->getFieldVal(pos);
  TMatrixD j(5, 5);
  // normal to plane
  TVector3 I = plane->getNormal();
  // U direction
  TVector3 J = plane->getU();
  // V direction
  TVector3 K = plane->getV();

  TVector3 H = B * (1. / B.Mag());
  TVector3 T = mom * (1. / mom.Mag());
  double alpha = (H.Cross(T)).Mag();

  TVector3 Z(0., 0., 1.);

  TVector3 U = Z.Cross(T) * (1. / (Z.Cross(T)).Mag());
  TVector3 V = T.Cross(U);
  TVector3 N = H.Cross(T) * (1. / alpha);
  double cosLambda = getCosLambda(mom);

  double Q = - B.Mag() * charge / mom.Mag();

  j(0, 0) = 1.;
  j(0, 1) = 0.;
  j(0, 2) = 0.;
  j(0, 3) = 0.;
  j(0, 4) = 0.;

  j(1, 0) = 0.;
  j(1, 1) = (T.Dot(I)) * (V.Dot(J));
  j(1, 2) = (T.Dot(I)) * (V.Dot(K));
  j(1, 3) = - alpha * Q * (T.Dot(J)) * (V.Dot(N));
  j(1, 4) = - alpha * Q * (T.Dot(K)) * (V.Dot(N));

  j(2, 0) = 0.;
  j(2, 1) = (T.Dot(I)) * (U.Dot(J)) / cosLambda;
  j(2, 2) = (T.Dot(I)) * (U.Dot(K)) / cosLambda;
  j(2, 3) = - alpha * Q * (T.Dot(J)) * (U.Dot(N)) / cosLambda;
  j(2, 4) = - alpha * Q * (T.Dot(K)) * (U.Dot(N)) / cosLambda;

  j(3, 0) = 0.;
  j(3, 1) = 0.;
  j(3, 2) = 0.;
  j(3, 3) = U.Dot(J);
  j(3, 4) = U.Dot(K);

  j(4, 0) = 0.;
  j(4, 1) = 0.;
  j(4, 2) = 0.;
  j(4, 3) = V.Dot(J);
  j(4, 4) = V.Dot(K);

  return j;
}

/**
 * @brief Projection matrix from Curvilinear to Local(measurement) coordinates
 * See paper from Strandlie & Wittek
 * Nuclear Instruments and Methods in Physics Research
 * Derivation of Jacobians for the propagation of covariance matrices of track parameters in homogeneous magnetic fields
 * A 566 (2006) 687-698
 * page 697, formula A.25 ... bottom right corner of the matrix
 * @param plane ...
 * @param mom ...
 * @return TMatrixD
 */
TMatrixD projCurv2Meas(SharedPlanePtr plane, TVector3 mom) {
  TMatrixD p(2, 2);

  //NOTE Following lines are version from Claus' python toy simulation
  /*TVector3 normal = plane->getNormal();
  double ey = normal[0];
  double ex = - normal[1];
  double ez = 1.;
  //proL2m[0][0] = HitHMatrix[0][3];
  //proL2m[1][1] = HitHMatrix[1][4];

  double cosLambda = getCosLambda(mom);
  double sinPhi = getSinPhi(mom);
  double sinLambda = sqrt(1. - cosLambda * cosLambda);
  double cosPhi = sqrt(1. - sinPhi * sinPhi);

  TMatrixD uvDir(2, 3);
  uvDir(0, 0) = - sinPhi; uvDir(0, 1) = cosPhi; uvDir(0, 2) = 0.;
  uvDir(1, 0) = -sinLambda * cosPhi; uvDir(1, 1) = -sinLambda * sinPhi; uvDir(1, 2) = cosLambda;

  TMatrixD mDir(2, 3);
  mDir(0, 0) = ex; mDir(0, 1) = ey; mDir(0, 2) = 0.;
  mDir(1, 0) = 0.; mDir(1, 1) = 0.; mDir(1, 2) = ez;

  p = uvDir * (mDir.T());
  return p.Invert();
  */

  // normal to plane
  TVector3 I = plane->getNormal();
  // U direction
  TVector3 J = plane->getU();
  // V direction
  TVector3 K = plane->getV();

  TVector3 T = mom * (1. / mom.Mag());
  TVector3 Z(0., 0., 1.);

  TVector3 U = Z.Cross(T) * (1. / (Z.Cross(T)).Mag());//
  TVector3 V = T.Cross(U);//

  double invTI = 1. / T.Dot(I); //

  p(0, 0) = V.Dot(K) * invTI;
  p(0, 1) = - U.Dot(K) * invTI;
  p(1, 0) = - V.Dot(J) * invTI;
  p(1, 1) = U.Dot(J) * invTI;

  return p;
}


void GFGbl::processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits) {
  //TGeoManager *gGeoManager;

  // flag for checking if desired sensor for alignment has been hitted
  bool flag = false;
  bool simple = false;

  TrackPoint* point_meas;
  AbsMeasurement* raw_meas;
  TVectorD measPrec(2);

  int npoints_meas = trk->getNumPointsWithMeasurement();

  std::vector<GblPoint> listOfPoints;
  unsigned int seedLabel = 0; // Add internal/external seed option in the future
  TMatrixDSym clCov(5), clSeed(5);
  TMatrixD jacPointToPoint(5, 5);

  jacPointToPoint.UnitMatrix();

  int n_gbl_points = 0;
  double Chi2 = 0.;
  int Ndf = 0;
  double lostWeight = 0.;

  genfit::FitStatus* fs = trk->getFitStatus(rep);
  genfit::KalmanFitStatus* kfs = dynamic_cast<genfit::KalmanFitStatus*>(fs);

  for (int ipoint_meas = 0; ipoint_meas < npoints_meas; ipoint_meas++) {
    point_meas = trk->getPointWithMeasurement(ipoint_meas);

    KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(point_meas->getFitterInfo(rep));
    // Current detector plane
    SharedPlanePtr plane = fi->getPlane();
    // StateOnPlane for Jacobian calculation and determination of length passed in sensor
    ReferenceStateOnPlane* reference = new ReferenceStateOnPlane(*fi->getReferenceState());
    //StateOnPlane* reference = new StateOnPlane(*fi->getFittedState());
    TVectorD state = reference->getState();

    // Get sensorId and material radiation length from geometry
    TVector3 o = plane->getO();
    gGeoManager->FindNode(o.x(), o.y(), o.z());
    int sensorId = getSensorIdFromPath(gGeoManager->GetPath());
    double radLength = gGeoManager->GetCurrentNode()->GetMedium()->GetMaterial()->GetRadLen();

    cout << "sensorId: " << sensorId << endl;
    raw_meas = point_meas->getRawMeasurement();
    TVectorD raw_coor = raw_meas->getRawHitCoords();

    // Here misalignment can be introduced to specific sensor
    if (sensorId == 8512) {
      raw_coor[0] += 0.100;
      raw_coor[1] -= 0.100;
    }

    // Covariance matrix of measurements
    const TMatrixDSym& raw_cov = raw_meas->getRawHitCov();
    // Projection matrix from repository state to measurement coords
    boost::scoped_ptr<const AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(rep));
    // Residual between measured position and reference track position
    TVectorD residual = raw_coor - HitHMatrix->Hv(state);

    // Extrapolate to middle of the sensor
    rep->extrapolateToPlane(*reference, plane, false, false);
    // I only get one half of the path inside sensor
    // ... should be good approximation
    double sensorHalfLenPassed = 0.;
    if (ipoint_meas < npoints_meas - 1) {
      // if not last point, extrapolate in track direction and stop at boundary
      sensorHalfLenPassed = rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas+1)->getFitterInfo(rep)->getPlane(), true, false);
    } else {
      // if last sensor, extrapolate to the previous
      // (and stop at boundary - of THIS sensor, not the one we extrapolate to!)
      sensorHalfLenPassed = rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas-1)->getFitterInfo(rep)->getPlane(), true, false);
    }
    // ensure positive length and double to get total path length inside sensor
    double sensorLenPassed = abs(sensorHalfLenPassed*2.);

    // For safety, extrapolate again to first sensor centre (no need to calculate Jacbian/Noise)
    // ... after this representation state is at this sensor plane
    rep->extrapolateToPlane(*reference, plane, false, false);

    // track direction at plane (in global coords)
    TVector3 trackDir = rep->getDir(*reference);
    // track momentum vector at plane (in global coords)
    TVector3 trackMom = rep->getMom(*reference);
    // track position vector at plane (in global coords)
    TVector3 trackPos = rep->getPos(*reference);
    // track momentum magnitude
    double trackMomMag = trackMom.Mag();
    // variable to store extrapolation length between this and next sensor
    double extLen = 1.;

    if (ipoint_meas < npoints_meas - 1) {
      TMatrixDSym noise;
      TVectorD deltaState;

      // Extrapolation to next plane with Jacobian calculation
      extLen = rep->extrapolateToPlane(*reference, trk->getPoint(ipoint_meas+1)->getFitterInfo(rep)->getPlane(), false, true);
      // Get the Jacobian from this plane to the following from the extrapolation
      rep->getForwardJacobianAndNoise(jacPointToPoint, noise, deltaState);
    }

    // Following code is to get jacobians for coordinate transformation
    /*
    // Jacobian local->curvilinear
    TMatrixD jacLoc2Cur(5, 5);
    // Jacobian curvilinear->local
    TMatrixD jacCur2Loc(5, 5);

    jacCur2Loc = jacLocal2Curvilinear(plane, trackMom, rep->getCharge(reference) / abs(rep->getCharge(reference)), trackPos);
    jacCur2Loc.Invert();
    if (ipoint_meas < npoints_meas -1){
    	// Now we are extrapolated at next plane, so get field and mom at that plane
    	// First get non-inverted jacobian of type local->curvilinear
    	jacLoc2Cur = jacLocal2Curvilinear(trk->getPoint(ipoint_meas+1)->getFitterInfo(rep)->getPlane(),
    																		 rep->getMom(reference),
    																     rep->getCharge(reference) / abs(rep->getCharge(reference)),
    																		 rep->getPos(reference)
    																		);
    }
    else
    {
    	jacCur2Loc.UnitMatrix();
    	jacLoc2Cur.UnitMatrix();
    }

    jacPointToPoint = jacLoc2Cur * (jacPointToPoint * jacCur2Loc);
    */

    // I now use only the simple jacobian
    if (simple) jacPointToPoint = getSimpleJac(extLen, 0.3*1.5/100, getCosLambda(trackMom));

    GblPoint point(jacPointToPoint);

    // Measurement precision = inverse covariance matrix (of sensor measurements)
    measPrec[0] = 1.0 / raw_cov(0,0);
    measPrec[1] = 1.0 / raw_cov(1,1);

    // Projection from local track coordinates (co-moving frame) to measurement coordinates
    TMatrixD proL2m(2, 2);

    // Projection from measurement to curvilinear track coordinates (co-moving frame)
    // ... old stuff

    if (simple) {
      // Get the projection from curvilinear sys. to meas. sys.
      proL2m = projCurv2Meas(plane, trackMom);
    }
    else {
      TMatrixD proM2L(2,2);
      proM2L.Zero();
      proM2L.Zero();
      proM2L(0,0) = HitHMatrix->getMatrix()(0,3);
      proM2L(1,1) = HitHMatrix->getMatrix()(1,4);
      proL2m = proM2L.Invert();
    }


    // Initial kinks = [0, 0]
    TVectorD scatResidual = TVectorD(2);
    scatResidual(0) = 0.;
    scatResidual(1) = 0.;
    // Passed fraction of radiation length in sensor
    double sensorXX0 = sensorLenPassed/radLength;

    //TVectorD scatPrecision = TVectorD(2);
    //scatPrecision(0) = 1./ getScattSigma2(trackMomMag, sensorXX0);
    //scatPrecision(1) = 1./ getScattSigma2(trackMomMag, sensorXX0) * getCosLambda2(trackMom);

    TMatrixDSym scatCov(2);
    if (simple) {
      scatCov(0, 0) = getScattSigma2(trackMomMag, sensorXX0);
      scatCov(1, 1) = getScattSigma2(trackMomMag, sensorXX0);
      scatCov(0, 1) = 0.;
      scatCov(1, 0) = 0.;
    }
    else {
      double c1 = trackDir.Dot(plane->getU());
      double c2 = trackDir.Dot(plane->getV());
      scatCov(0, 0) = 1. - c2*c2;
      scatCov(1, 1) = 1. - c1*c1;
      scatCov(0, 1) = c1*c2;
      scatCov(1, 0) = c1*c2;
      scatCov *= getScattSigma2(trackMomMag, sensorXX0) / (1. - c1 * c1 - c2 * c2) / (1. - c1 * c1 - c2 * c2) ;
    }
    point.addMeasurement(proL2m, residual, measPrec);
    point.addScatterer(scatResidual, scatCov.Invert());
    //point.addScatterer(scatResidual, scatPrecision);

    //Add global derivatives to the point

    // sensor label = sensorID * 10, then last digit is label for global derivative for the sensor
    int label = sensorId * 10;
    cout << label << endl;
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

    //WARNING This now aligns olny one specified sensor
    if (sensorId == 8512 && kfs->getBackwardNdf() == 7)	{
      flag = true;
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

      //TODO Global derivatives continuation
      // ... needs changes in labels and positions of derivatives values in matrix
      /*
      // Global derivatives for movement of whole detector system in global coordinates 
      //TODO: Usage of this would require Hierarchy Constraints to be provided to MP2

      // sensor centre position in global system
      TVector3 detPos = plane->getO();
      cout << "detPos" << endl;
      detPos.Print();

      // global prediction from raw measurement
      TVector3 pred = detPos + uPos * uDir + vPos * vDir;
      cout << "pred" << endl;
      pred.Print();

      double xPred = pred[0];
      double yPred = pred[1];
      double zPred = pred[2];

      // scalar product of sensor normal and track direction
      double tn = tDir.Dot(nDir);
      cout << "tn" << endl;
      cout << tn << endl;

      // derivatives of local residuals versus measurements
      TMatrixD drdm(3, 3);
      drdm.UnitMatrix();
      for(int row = 0; row < 3; row++)
      	for(int col = 0; col < 3; col++)
      		drdm(row, col) -= tDir[row] * nDir[col] / tn;
      	
      cout << "drdm" << endl;
      drdm.Print();

      // derivatives of measurements versus global alignment parameters
      TMatrixD dmdg(3, 6);
      dmdg.Zero();
      dmdg(0, 0) = 1.; dmdg(0, 4) = -zPred; dmdg(0, 5) = yPred;
      dmdg(1, 1) = 1.; dmdg(1, 3) = zPred;  dmdg(1, 5) = -xPred;
      dmdg(2, 2) = 1.; dmdg(2, 3) = -yPred; dmdg(2, 4) = xPred;

      cout << "dmdg" << endl;
      dmdg.Print();

      // derivatives of local residuals versus global alignment parameters
      TMatrixD drldrg(3, 3);
      drldrg.Zero();
      drldrg(0, 0) = uDir[0]; drldrg(0, 1) = uDir[1]; drldrg(0, 2) = uDir[2];
      drldrg(1, 0) = vDir[0]; drldrg(1, 1) = vDir[1]; drldrg(1, 2) = vDir[2];

      cout << "drldrg" << endl;
      drldrg.Print();

      cout << "drdm * dmdg" << endl;
      (drdm * dmdg).Print();

      // derivatives of local residuals versus rigid body parameters
      TMatrixD drldg(3, 6);
      drldg = drldrg * (drdm * dmdg);

      cout << "drldg" << endl;
      drldg.Print();

      // offset to determine labels for sensor sets or individual layers
      // 0: PXD, TODO 1: SVD, or individual layers
      int offset = 0;

      derGlobal(0, 5) = drldg(0, 0); labGlobal(0, 5) = offset + 1;
      derGlobal(0, 6) = drldg(0, 1); labGlobal(0, 6) = offset + 2;
      derGlobal(0, 7) = drldg(0, 2); labGlobal(0, 7) = offset + 3;
      derGlobal(0, 8) = drldg(0, 3);	labGlobal(0, 8) = offset + 4;
      derGlobal(0, 9) = drldg(0, 4); labGlobal(0, 9) = offset + 5;
      derGlobal(0,10) = drldg(0, 5);	labGlobal(0,10) = offset + 6;		
      	
      derGlobal(1, 5) = drldg(1, 0); labGlobal(1, 5) = offset + 1;
      derGlobal(1, 6) = drldg(1, 1); labGlobal(1, 6) = offset + 2;
      derGlobal(1, 7) = drldg(1, 2); labGlobal(1, 7) = offset + 3;
      derGlobal(1, 8) = drldg(1, 3);	labGlobal(1, 8) = offset + 4;
      derGlobal(1, 9) = drldg(1, 4); labGlobal(1, 9) = offset + 5;
      derGlobal(1,10) = drldg(1, 5);	labGlobal(1,10) = offset + 6;			
      */

      point.addGlobals(labGlobal, derGlobal);
    }

    listOfPoints.push_back(point);
    n_gbl_points++;
    // Free memory on the heap
    delete reference;
  }
  if (flag && kfs->getBackwardNdf() == 7) {
    //if (n_gbl_points >= 2 && flag) {
    GblTrajectory * traj = 0;
    try {
      traj = new GblTrajectory(listOfPoints, seedLabel, clSeed);
      traj->fit(Chi2, Ndf, lostWeight);
    } catch(...) {
      // Gbl failed critically (usually GblPoint::getDerivatives ... singular matrix inversion)
      delete traj;
      return;
    }
    /*
    //TODO Second GBL external iteration ... for electrons with Brehmstrahlung
    for (unsigned int p = 0; p < listOfPoints.size(); p++){
    	unsigned int label = p + 1;
    	unsigned int numRes;
    	TVectorD residuals(2);
    	TVectorD measErrors(2);
    	TVectorD resErrors(2);
    	TVectorD downWeights(2);

    	traj->getScatResults(label, numRes, residuals, measErrors, resErrors, downWeights);
    	resErrors[0] = 1. / resErrors[0] / resErrors[0];
    	resErrors[1] = 1. / resErrors[1] / resErrors[1];
    	// TODO: in case of more iterations here must be also kinks already added to points
    	// in previous iteration
    	residuals[0] = - residuals[0];
    	residuals[1] = - residuals[1];
    	listOfPoints[p].addScatterer(residuals, resErrors);
    }

    delete traj;
    Chi2 = 0.; Ndf = 0; lostWeight = 0.;

    try{
    	traj = new GblTrajectory(listOfPoints, seedLabel, clSeed);
    	traj->fit(Chi2, Ndf, lostWeight, "C");
    }catch(...){
    	// Gbl failed critically (usually GblPoint::getDerivatives ... singular matrix inversion
    	delete traj;
    	return;
    }*/

    traj->milleOut(*milleFile);
    delete traj;

#ifdef DEBUG
    ofs_output << kfs->getBackwardNdf() << " " << kfs->getBackwardChi2() << " "
	       << Ndf << " " << Chi2 << std::endl;
#endif
  }

}
