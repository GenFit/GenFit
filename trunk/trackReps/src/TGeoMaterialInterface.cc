/* Copyright 2008-2014, Technische Universitaet Muenchen,
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

#include <TGeoMaterialInterface.h>
#include <Exception.h>

#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoManager.h>
#include <assert.h>
#include <math.h>

static const bool debug = false;
//static const bool debug = true;

namespace genfit {

double MeanExcEnergy_get(int Z);
double MeanExcEnergy_get(TGeoMaterial*);


bool
TGeoMaterialInterface::initTrack(double posX, double posY, double posZ,
                                   double dirX, double dirY, double dirZ){
  #ifdef DEBUG
  std::cout << "TGeoMaterialInterface::initTrack. \n";
  std::cout << "Pos    "; TVector3(posX, posY, posZ).Print();
  std::cout << "Dir    "; TVector3(dirX, dirY, dirZ).Print();
  #endif

  // Move to the new point.
  bool result = !gGeoManager->IsSameLocation(posX, posY, posZ, kTRUE);
  // Set the intended direction.
  gGeoManager->SetCurrentDirection(dirX, dirY, dirZ);
  return result;
}


void
TGeoMaterialInterface::getMaterialParameters(double& density,
                                               double& Z,
                                               double& A,
                                               double& radiationLength,
                                               double& mEE){

  TGeoMaterial* mat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

  density         = mat->GetDensity();
  Z               = mat->GetZ();
  A               = mat->GetA();
  radiationLength = mat->GetRadLen();
  mEE             = MeanExcEnergy_get(mat);

}


void
TGeoMaterialInterface::getMaterialParameters(MaterialProperties& parameters) {

  TGeoMaterial* mat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();

  parameters.setMaterialProperties(mat->GetDensity(),
      mat->GetZ(),
      mat->GetA(),
      mat->GetRadLen(),
      MeanExcEnergy_get(mat));

}


double
TGeoMaterialInterface::findNextBoundary(const RKTrackRep* rep,
                                          const M1x7& stateOrig,
                                          double sMax, // signed
                                          bool varField)
{
  const double delta(1.E-2); // cm, distance limit beneath which straight-line steps are taken.
  const double epsilon(1.E-1); // cm, allowed upper bound on arch
			       // deviation from straight line

  M1x3 SA;
  M1x7 state7, oldState7;
  memcpy(oldState7, stateOrig, sizeof(state7));

  int stepSign(sMax < 0 ? -1 : 1);

  double s = 0;  // trajectory length to boundary

  const unsigned maxIt = 300;
  unsigned it = 0;

  // Initialize the geometry to the current location (set by caller).
  gGeoManager->FindNextBoundary(fabs(sMax) - s);
  double safety = gGeoManager->GetSafeDistance();
  double slDist = gGeoManager->GetStep();
  double step = slDist;

  while (1) {
    if (++it > maxIt){
      Exception exc("TGeoMaterialInterface::findNextBoundary ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // No boundary in sight?
    if (s + safety > fabs(sMax)) {
      if (debug)
	std::cout << "   next boundary is further away than sMax \n";
      return stepSign*(s + safety); //sMax;
    }

    // Are we at the boundary?
    if (slDist < delta) {
      if (debug)
	std::cout << "   very close to the boundary -> return"
		  << " stepSign*(s + slDist) = "
		  << stepSign << "*(" << s + slDist << ")\n";
      return stepSign*(s + slDist);
    }

    // We have to find whether there's any boundary on our path.

    // Follow curved arch, then see if we may have missed a boundary.
    // Always propagate complete way from original start to avoid
    // inconsistent extrapolations.
    memcpy(state7, stateOrig, sizeof(state7));
    rep->RKPropagate(state7, NULL, SA, stepSign*(s + step), varField);

    // Straight line distance² between extrapolation finish and
    // the end of the previously determined safe segment.
    double dist2 = (pow(state7[0] - oldState7[0], 2)
		    + pow(state7[1] - oldState7[1], 2)
		    + pow(state7[2] - oldState7[2], 2));
    // Maximal lateral deviation².
    double maxDeviation2 = 0.25*(step*step - dist2);

    if (step > safety
	&& maxDeviation2 > epsilon*epsilon) {
      // Need to take a shorter step to reliably estimate material,
      // but only if we didn't move by safety.  In order to avoid
      // issues with roundoff we check 'step > safety' instead of
      // 'step != safety'.  If we moved by safety, there couldn't have
      // been a boundary that we missed along the path, but we could
      // be on the boundary.

      // Take a shorter step, but never shorter than safety.
      step = std::max(step / 2, safety);
    } else {
      gGeoManager->PushPoint();
      bool volChanged = initTrack(state7[0], state7[1], state7[2],
				  stepSign*state7[3], stepSign*state7[4],
				  stepSign*state7[5]);

      if (volChanged) {
	// Move back to start.
	gGeoManager->PopPoint();

	// Extrapolation may not take the exact step length we asked
	// for, so it can happen that a requested step < safety takes
	// us across the boundary.  This is then the best estimate we
	// can get of the distance to the boundary with the stepper.
	if (step <= safety)
	  return stepSign*(s + step);

	// Volume changed during the extrapolation.  Take a shorter
	// step, but never shorter than safety.
	step = std::max(step / 2, safety);
      } else {
	// we're in the new place, the step was safe, advance
	s += step;

	memcpy(oldState7, state7, sizeof(state7));
	gGeoManager->PopDummy();  // Pop stack, but stay in place.

	gGeoManager->FindNextBoundary(fabs(sMax) - s);
	safety = gGeoManager->GetSafeDistance();
	step = slDist = gGeoManager->GetStep();
      }
    }
  }
}


double
TGeoMaterialInterface::findNextBoundaryAndStepStraight(double sMax) {

  gGeoManager->FindNextBoundaryAndStep(sMax);
  return gGeoManager->GetStep();

}




/*
Reference for elemental mean excitation energies at:
http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html

Code ported from GEANT 3
*/

const int MeanExcEnergy_NELEMENTS = 93; // 0 = vacuum, 1 = hydrogen, 92 = uranium
const double MeanExcEnergy_vals[] = {1.e15, 19.2, 41.8, 40.0, 63.7, 76.0, 78., 82.0, 95.0, 115.0, 137.0, 149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0};


double
MeanExcEnergy_get(int Z) {
  assert(Z >= 0 && Z < MeanExcEnergy_NELEMENTS);
  return MeanExcEnergy_vals[Z];
}


double
MeanExcEnergy_get(TGeoMaterial* mat) {
  if (mat->IsMixture()) {
    double logMEE = 0.;
    double denom  = 0.;
    TGeoMixture* mix = (TGeoMixture*)mat;
    for (int i = 0; i < mix->GetNelements(); ++i) {
      int index = int(floor((mix->GetZmixt())[i]));
      logMEE += 1. / (mix->GetAmixt())[i] * (mix->GetWmixt())[i] * (mix->GetZmixt())[i] * log(MeanExcEnergy_get(index));
      denom  += (mix->GetWmixt())[i] * (mix->GetZmixt())[i] * 1. / (mix->GetAmixt())[i];
    }
    logMEE /= denom;
    return exp(logMEE);
  } else { // not a mixture
    int index = int(floor(mat->GetZ()));
    return MeanExcEnergy_get(index);
  }
}


} /* End of namespace genfit */
