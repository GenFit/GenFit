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

#include "MeasurementCreator.h"

#include <iostream>

#include <PlanarMeasurement.h>
#include <ProlateSpacepointMeasurement.h>
#include <SpacepointMeasurement.h>
#include <WireMeasurement.h>
#include <WirePointMeasurement.h>

#include <TRandom.h>
#include <TMath.h>

#include <assert.h>
#include <math.h>

namespace genfit {

MeasurementCreator::MeasurementCreator() :
    trackModel_(NULL),
    resolution_(0.01),
    resolutionWire_(0.1),
    outlierProb_(0),
    outlierRange_(2),
    thetaDetPlane_(90),
    phiDetPlane_(0),
    wireCounter_(0),
    wireDir_(0.,0.,1.),
    minDrift_(0),
    maxDrift_(2),
    idealLRResolution_(true),
    useSkew_(false),
    skewAngle_(5),
    nSuperLayer_(5),
    measurementCounter_(0),
    debug_(false)
{
  ;
}


std::vector<genfit::AbsMeasurement*> MeasurementCreator::create(eMeasurementType type, double tracklength, bool& outlier, int& lr) {

  outlier = false;
  lr = 0;
  std::vector<AbsMeasurement*> retVal;
  genfit::AbsMeasurement* measurement;

  TVector3 point, dir;
  trackModel_->getPosDir(tracklength, point, dir);


  TVector3 planeNorm(dir);
  planeNorm.SetTheta(thetaDetPlane_*TMath::Pi()/180);
  planeNorm.SetPhi(planeNorm.Phi()+phiDetPlane_);
  static const TVector3 z(0,0,1);
  static const TVector3 x(1,0,0);


  TVector3 currentWireDir(wireDir_);
  TVector3 wirePerp;

  if (type == Wire ||
      type == WirePoint){

    // skew layers
    if (useSkew_ && (int)((double)wireCounter_/(double)nSuperLayer_)%2 == 1) {
      TVector3 perp(wireDir_.Cross(dir));
      if ((int)((double)wireCounter_/(double)nSuperLayer_)%4 == 1){
        currentWireDir.Rotate(skewAngle_*TMath::Pi()/180, wireDir_.Cross(perp));
      }
      else currentWireDir.Rotate(-skewAngle_*TMath::Pi()/180, wireDir_.Cross(perp));
    }
    currentWireDir.SetMag(1.);

    // left/right
    lr = 1;
    wirePerp = dir.Cross(currentWireDir);
    if (gRandom->Uniform(-1,1) >= 0) {
      wirePerp *= -1.;
      lr = -1;
    }
    wirePerp.SetMag(gRandom->Uniform(minDrift_, maxDrift_));
  }

  if (outlierProb_ > gRandom->Uniform(1.)) {
    outlier = true;
    if(debug_)  std::cerr << "create outlier" << std::endl;
  }


  switch(type){
    case Pixel: {
      if (debug_) std::cerr << "create PixHit" << std::endl;

      genfit::SharedPlanePtr plane(new genfit::DetPlane(point, planeNorm.Cross(z), (planeNorm.Cross(z)).Cross(planeNorm)));

      TVectorD hitCoords(2);
      if (outlier) {
        hitCoords(0) = gRandom->Uniform(-outlierRange_, outlierRange_);
        hitCoords(1) = gRandom->Uniform(-outlierRange_, outlierRange_);
      }
      else {
        hitCoords(0) = gRandom->Gaus(0,resolution_);
        hitCoords(1) = gRandom->Gaus(0,resolution_);
      }

      TMatrixDSym hitCov(2);
      hitCov(0,0) = resolution_*resolution_;
      hitCov(1,1) = resolution_*resolution_;

      measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, measurementCounter_);
      retVal.push_back(measurement);
    }
    break;

    case Spacepoint: {
      if (debug_) std::cerr << "create SpacepointHit" << std::endl;

      TVectorD hitCoords(3);
      if (outlier) {
        hitCoords(0) = gRandom->Uniform(point.X()-outlierRange_, point.X()+outlierRange_);
        hitCoords(1) = gRandom->Uniform(point.Y()-outlierRange_, point.Y()+outlierRange_);
        hitCoords(2) = gRandom->Uniform(point.Z()-outlierRange_, point.Z()+outlierRange_);
      }
      else {
        hitCoords(0) = gRandom->Gaus(point.X(),resolution_);
        hitCoords(1) = gRandom->Gaus(point.Y(),resolution_);
        hitCoords(2) = gRandom->Gaus(point.Z(),resolution_);
      }

      TMatrixDSym hitCov(3);
      hitCov(0,0) = resolution_*resolution_;
      hitCov(1,1) = resolution_*resolution_;
      hitCov(2,2) = resolution_*resolution_;

      measurement = new genfit::SpacepointMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      retVal.push_back(measurement);
    }
    break;

    case ProlateSpacepoint: {
      if (debug_) std::cerr << "create ProlateSpacepointHit" << std::endl;

      TVectorD hitCoords(3);
      if (outlier) {
        hitCoords(0) = gRandom->Uniform(point.X()-outlierRange_, point.X()+outlierRange_);
        hitCoords(1) = gRandom->Uniform(point.Y()-outlierRange_, point.Y()+outlierRange_);
        hitCoords(2) = gRandom->Uniform(point.Z()-outlierRange_, point.Z()+outlierRange_);
      }
      else {
        hitCoords(0) = point.X();
        hitCoords(1) = point.Y();
        hitCoords(2) = point.Z();
      }

      TMatrixDSym hitCov(3);
      hitCov(0,0) = resolution_*resolution_;
      hitCov(1,1) = resolution_*resolution_;
      hitCov(2,2) = resolutionWire_*resolutionWire_;

      // rotation matrix
      TVector3 xp = currentWireDir.Orthogonal();
      xp.SetMag(1);
      TVector3 yp = currentWireDir.Cross(xp);
      yp.SetMag(1);

      TMatrixD rot(3,3);

      rot(0,0) = xp.X();  rot(0,1) = yp.X();  rot(0,2) = currentWireDir.X();
      rot(1,0) = xp.Y();  rot(1,1) = yp.Y();  rot(1,2) = currentWireDir.Y();
      rot(2,0) = xp.Z();  rot(2,1) = yp.Z();  rot(2,2) = currentWireDir.Z();

      // smear
      TVectorD smearVec(3);
      smearVec(0) = resolution_;
      smearVec(1) = resolution_;
      smearVec(2) = resolutionWire_;
      smearVec *= rot;
      if (!outlier) {
        hitCoords(0) += gRandom->Gaus(0, smearVec(0));
        hitCoords(1) += gRandom->Gaus(0, smearVec(1));
        hitCoords(2) += gRandom->Gaus(0, smearVec(2));
      }


      // rotate cov
      hitCov.Similarity(rot);

      measurement = new genfit::ProlateSpacepointMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      static_cast<genfit::ProlateSpacepointMeasurement*>(measurement)->setLargestErrorDirection(currentWireDir);
      retVal.push_back(measurement);
    }
    break;

    case StripU: case StripV: case StripUV : {
      if (debug_) std::cerr << "create StripHit" << std::endl;

      TVector3 vU, vV;
      vU = planeNorm.Cross(z);
      vV = (planeNorm.Cross(z)).Cross(planeNorm);
      genfit::SharedPlanePtr plane(new genfit::DetPlane(point, vU, vV));

      TVectorD hitCoords(1);
      if (outlier)
        hitCoords(0) = gRandom->Uniform(-outlierRange_, outlierRange_);
      else
        hitCoords(0) = gRandom->Gaus(0,resolution_);

      TMatrixDSym hitCov(1);
      hitCov(0,0) = resolution_*resolution_;

      measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, measurementCounter_);
      if (type == StripV)
        static_cast<genfit::PlanarMeasurement*>(measurement)->setStripV();
      retVal.push_back(measurement);


      if (type == StripUV) {
        if (outlier)
          hitCoords(0) = gRandom->Uniform(-outlierRange_, outlierRange_);
        else
          hitCoords(0) = gRandom->Gaus(0,resolution_);

        hitCov(0,0) = resolution_*resolution_;

        measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
        static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, measurementCounter_);
        static_cast<genfit::PlanarMeasurement*>(measurement)->setStripV();
        retVal.push_back(measurement);
      }
    }
    break;

    case Wire: {
      if (debug_) std::cerr << "create WireHit" << std::endl;

      if (outlier) {
        wirePerp.SetMag(gRandom->Uniform(outlierRange_));
      }

      TVectorD hitCoords(7);
      hitCoords(0) = (point-wirePerp-currentWireDir).X();
      hitCoords(1) = (point-wirePerp-currentWireDir).Y();
      hitCoords(2) = (point-wirePerp-currentWireDir).Z();

      hitCoords(3) = (point-wirePerp+currentWireDir).X();
      hitCoords(4) = (point-wirePerp+currentWireDir).Y();
      hitCoords(5) = (point-wirePerp+currentWireDir).Z();

      if (outlier)
        hitCoords(6) = gRandom->Uniform(outlierRange_);
      else
        hitCoords(6) = gRandom->Gaus(wirePerp.Mag(), resolution_);

      TMatrixDSym hitCov(7);
      hitCov(6,6) = resolution_*resolution_;


      measurement = new genfit::WireMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      if (idealLRResolution_){
        static_cast<genfit::WireMeasurement*>(measurement)->setLeftRightResolution(lr);
      }
      ++wireCounter_;
      retVal.push_back(measurement);
    }
    break;

    case WirePoint: {
      if (debug_) std::cerr << "create WirePointHit" << std::endl;

      if (outlier) {
        wirePerp.SetMag(gRandom->Uniform(outlierRange_));
      }

      TVectorD hitCoords(8);
      hitCoords(0) = (point-wirePerp-currentWireDir).X();
      hitCoords(1) = (point-wirePerp-currentWireDir).Y();
      hitCoords(2) = (point-wirePerp-currentWireDir).Z();

      hitCoords(3) = (point-wirePerp+currentWireDir).X();
      hitCoords(4) = (point-wirePerp+currentWireDir).Y();
      hitCoords(5) = (point-wirePerp+currentWireDir).Z();

      if (outlier) {
        hitCoords(6) = gRandom->Uniform(outlierRange_);
        hitCoords(7) = gRandom->Uniform(currentWireDir.Mag()-outlierRange_, currentWireDir.Mag()+outlierRange_);
      }
      else {
        hitCoords(6) = gRandom->Gaus(wirePerp.Mag(), resolution_);
        hitCoords(7) = gRandom->Gaus(currentWireDir.Mag(), resolutionWire_);
      }


      TMatrixDSym hitCov(8);
      hitCov(6,6) = resolution_*resolution_;
      hitCov(7,7) = resolutionWire_*resolutionWire_;

      measurement = new genfit::WirePointMeasurement(hitCoords, hitCov, int(type), measurementCounter_, nullptr);
      if (idealLRResolution_){
        static_cast<genfit::WirePointMeasurement*>(measurement)->setLeftRightResolution(lr);
      }
      ++wireCounter_;
      retVal.push_back(measurement);
    }
    break;

    default:
      std::cerr << "measurement type not defined!" << std::endl;
      exit(0);
  }

  return retVal;

}


void MeasurementCreator::reset() {
  wireCounter_ = 0;
  measurementCounter_ = 0;
}

} /* End of namespace genfit */
