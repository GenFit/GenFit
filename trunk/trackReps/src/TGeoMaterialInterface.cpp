/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

//#define DEBUG


namespace genfit {

double MeanExcEnergy_get(int Z);
double MeanExcEnergy_get(TGeoMaterial*);


TGeoMaterialInterface::TGeoMaterialInterface() : whichNavig_(0) {
  gGeoManager->AddNavigator();
}

void
TGeoMaterialInterface::initTrack(double posX, double posY, double posZ,
                                   double dirX, double dirY, double dirZ){
  #ifdef DEBUG
  std::cout << "TGeoMaterialInterface::initTrack. \n";
  std::cout << "Pos    "; TVector3(posX, posY, posZ).Print();
  std::cout << "Dir    "; TVector3(dirX, dirY, dirZ).Print();
  #endif

  const Double_t * pt = gGeoManager->GetCurrentPoint();

  if (posX == pt[0] &&
      posY == pt[1] &&
      posZ == pt[2]) {
    // position does not change
    gGeoManager->SetCurrentDirection(dirX, dirY, dirZ);
    #ifdef DEBUG
    std::cout << "just init dir (default navigator)! \n";
    #endif
  }
  else {

    whichNavig_ = (whichNavig_ + 1) % 2;
    gGeoManager->SetCurrentNavigator(whichNavig_);
    pt = gGeoManager->GetCurrentPoint();

    if (posX == pt[0] &&
        posY == pt[1] &&
        posZ == pt[2]) {
      // position does not change
      gGeoManager->SetCurrentDirection(dirX, dirY, dirZ);
      #ifdef DEBUG
      std::cout << "just init dir (second navigator)! \n";
      #endif
    }
    else {
      gGeoManager->InitTrack(posX, posY, posZ,
                             dirX, dirY, dirZ);
    }
  }

  #ifdef DEBUG
  std::cout << "\n";
  #endif

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
                                          bool varField){

  const double delta(1.E-2); // cm
  double s(0), safety(0), slDist(0);
  M1x3 SA;
  M1x7 state7;
  memcpy(state7, stateOrig, sizeof(stateOrig));

  int stepSign(1);
  if (sMax < 0) stepSign = -1;

  unsigned int maxIt(300), it(0);


  while (true) {

    if (++it > maxIt){
      Exception exc("TGeoMaterialInterface::findNextBoundary ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    safety = gGeoManager->Safety(); // unsigned; distance to closest boundary

#ifdef DEBUG
    std::cout << "   TGeoMaterialInterface::findNextBoundary: Iteration " << it << ". Safety = " << safety << ". slDist = " << slDist << ". Step so far = " << s << "\n";
    std::cout << "   Material before step: " << gGeoManager->GetCurrentVolume()->GetMedium()->GetName() << "\n";
#endif

    if (fabs(s + stepSign*safety) > fabs(sMax)) { // next boundary is further away than sMax
#ifdef DEBUG
      std::cout << "   next boundary is further away than sMax \n";
#endif
      return s + stepSign*safety;
    }

    gGeoManager->FindNextBoundary();
    slDist = gGeoManager->GetStep(); // unsigned; straight line distance to next boundary along step direction

    if (slDist < delta) { // very near the boundary
#ifdef DEBUG
      std::cout << "   very near the boundary -> return s + stepSign*slDist; = " << s + stepSign*slDist << "\n";
#endif
      return s + stepSign*slDist;
    }
    else {
#ifdef DEBUG
      std::cout << "   make RKutta step \n";
#endif
      // check if we would cross a boundary when making slDist step
      double tryStep = 0.9 * stepSign*slDist;

      bool safe = false;
      if (fabs(tryStep) < fabs(safety))
        safe = true;

      if (!safe) {
        memcpy(state7, stateOrig, sizeof(state7)); // propagate complete way from original start
        rep->RKPropagate(state7, NULL, SA, s+tryStep, varField);
        // init for checking
        initTrack(stateOrig[0], stateOrig[1], stateOrig[2],
                  state7[0]-stateOrig[0], state7[1]-stateOrig[1], state7[2]-stateOrig[2]);
      }
      // check: gGeoManager->GetStep() gives us the max sl step size from stateOrig to state7
      if (!safe && gGeoManager->GetStep() > fabs(tryStep)) {
        s += tryStep;
        // init for next iteration
        initTrack(state7[0], state7[1], state7[2],  stepSign*state7[3], stepSign*state7[4], stepSign*state7[5]);
        #ifdef DEBUG
          std::cout << "   tried and its safe to make a step of  " << stepSign*tryStep << "\n";
        #endif
      }
      else { // step along safety
        s += stepSign*safety;
        memcpy(state7, stateOrig, sizeof(state7)); // propagate complete way from original start
        rep->RKPropagate(state7, NULL, SA, s, varField);
        // init for next iteration
        initTrack(state7[0], state7[1], state7[2],  stepSign*state7[3], stepSign*state7[4], stepSign*state7[5]);
#ifdef DEBUG
    std::cout << "   step along safety  " << stepSign*safety << "\n";
#endif
      }

    }

#ifdef DEBUG
    std::cout << "   Material after step: " << gGeoManager->GetCurrentVolume()->GetMedium()->GetName() << "\n";
#endif

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
