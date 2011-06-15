/* Copyright 2008-2009, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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

#include"GFMaterialEffects.h"
#include<iostream>
#include "stdlib.h"

#include"TDatabasePDG.h"
#include"TGeoMaterial.h"
#include"TGeoManager.h"

#include"math.h"
#include"assert.h"


GFMaterialEffects* GFMaterialEffects::finstance = NULL;

float MeanExcEnergy_get(int Z);
float MeanExcEnergy_get(TGeoMaterial*);


GFMaterialEffects::GFMaterialEffects():
  fEnergyLossBetheBloch(true), fNoiseBetheBloch(true),
  fNoiseCoulomb(true),
  fEnergyLossBrems(true), fNoiseBrems(true),
  me(0.510998910E-3){
}

GFMaterialEffects::~GFMaterialEffects(){
}

GFMaterialEffects* GFMaterialEffects::getInstance() {
  if(finstance == NULL) finstance = new GFMaterialEffects();
  return finstance;
}

void GFMaterialEffects::destruct() {
  if(finstance != NULL) {
    delete finstance;
    finstance = NULL;
  }
}


double GFMaterialEffects::effects(const std::vector<TVector3>& points, 
                                  const std::vector<double>& pointPaths, 
                                  const double& mom,
                                  const int& pdg,
                                  const bool& doNoise,
                                        TMatrixT<double>* noise,
                                  const TMatrixT<double>* jacobian,
                                  const TVector3* directionBefore, 
                                  const TVector3* directionAfter){

  assert(points.size()==pointPaths.size());

  double momLoss=0.;
                     
  for(unsigned int i=1;i<points.size();++i){
    TVector3 p1=points.at(i-1);
    TVector3 p2=points.at(i);
    TVector3 dir=p2-p1;
    double dist=dir.Mag();
    double realPath = pointPaths.at(i);
    
    if (dist > 1.E-8) { // do material effects only if distance is not too small
      dir*=1./dist; //normalize dir

      double X(0.);
      double matDensity, matZ, matA, radiationLength, mEE;
      double step;
      
      gGeoManager->InitTrack(p1.X(),p1.Y(),p1.Z(),dir.X(),dir.Y(),dir.Z());

      while(X<dist){
        
        gGeoManager->FindNextBoundaryAndStep(dist-X);
        step = gGeoManager->GetStep();
        
        assert(gGeoManager->GetCurrentVolume()->GetMedium()!=NULL);
        TGeoMaterial * mat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();
        
        matZ = mat->GetZ();
        
        if(matZ>1.E-3){ // don't calculate energy loss for vacuum        

          matDensity      = mat->GetDensity();
          matA            = mat->GetA();
          radiationLength = mat->GetRadLen();
          mEE             = MeanExcEnergy_get(mat);

          if (fEnergyLossBetheBloch)
            momLoss += realPath/dist * this->energyLossBetheBloch(step, mom, pdg, matDensity, matZ, matA, mEE);
          if (fEnergyLossBetheBloch && fNoiseBetheBloch)
            this->noiseBetheBloch(step, mom, pdg, matDensity, matZ, matA, mEE, noise);

          if (fNoiseCoulomb)
            this->noiseCoulomb(step, mom, pdg, matZ, radiationLength, noise, jacobian, directionBefore, directionAfter);

          if (fEnergyLossBrems)
            momLoss += realPath/dist * this->energyLossBrems(step, mom, pdg, matDensity, matZ, matA, radiationLength);
          if (fEnergyLossBrems && fNoiseBrems)
            this->noiseBrems(step, mom, pdg, radiationLength, noise);
        }
        X += step;
      }
    }
  }

  return momLoss;
}


double GFMaterialEffects::stepper(const double& maxDist,
                                  const double& posx,
                                  const double& posy,
                                  const double& posz,
                                  const double& dirx,
                                  const double& diry,
                                  const double& dirz,
                                  const double& mom,
                                  const int& pdg){

  static const double maxPloss = .005; // maximum relative momentum loss allowed

  gGeoManager->InitTrack(posx,posy,posz,dirx,diry,dirz);

  double X(0.);
  double dP = 0.;
  double momLoss = 0.;
  double matDensity, matZ, matA, radiationLength, mEE;
  double step;

  while(X<maxDist){
    
    gGeoManager->FindNextBoundaryAndStep(maxDist-X);
    step = gGeoManager->GetStep();
    
    assert(gGeoManager->GetCurrentVolume()->GetMedium()!=NULL);
    TGeoMaterial * mat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();
    
    matZ = mat->GetZ();
      
    if(matZ>1.E-3){ // don't calculate energy loss for vacuum       
 
      matDensity      = mat->GetDensity();
      matA            = mat->GetA();
      radiationLength = mat->GetRadLen();
      mEE             = MeanExcEnergy_get(mat);

      if (fEnergyLossBetheBloch)
        momLoss += this->energyLossBetheBloch(step, mom, pdg, matDensity, matZ, matA, mEE);

      if (fEnergyLossBrems)
        momLoss += this->energyLossBrems(step, mom, pdg, matDensity, matZ, matA, radiationLength);
    }
    
    if(dP + momLoss > mom*maxPloss){
      double fraction = (mom*maxPloss-dP)/momLoss;
      assert(fraction>0.&&fraction<1.);
      dP+=fraction*momLoss;
      X+=fraction*step;
      break;
    }
    
    dP += momLoss;
    X += step;
  }

  return X;                 
}


void GFMaterialEffects::getParticleParameters(const int&    pdg,
                                            double& charge,
                                            double& mass){
  TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(pdg);
  charge = part->Charge()/(3.);
  mass = part->Mass();
};

double GFMaterialEffects::getParticleMass (const int& pdg){
  TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(pdg);
  return part->Mass();
}



//---- Energy-loss and Noise calculations -----------------------------------------

double GFMaterialEffects::energyLossBetheBloch(const double& step,
                                               const double& mom,
                                               const int&    pdg,
                                               const double& matDensity,
                                               const double& matZ,
                                               const double& matA,
                                               const double& meanExcitationEnergy){

  double charge, mass;
  getParticleParameters(pdg, charge, mass);

  const double beta = mom/sqrt(mass*mass+mom*mom);
  double dedx = 0.307075*matZ/matA*matDensity/(beta*beta)*charge*charge;

  //for numerical stability
  double gammaSquare = 1.-beta*beta;
  if(gammaSquare>1.E-10) gammaSquare = 1./gammaSquare;
  else gammaSquare = 1.E10;
  double gamma = sqrt(gammaSquare);

  double massRatio = me/mass;
  double argument = gammaSquare*beta*beta*me*1.E3*2./((1.E-6*meanExcitationEnergy) * sqrt(1+2*sqrt(gammaSquare)*massRatio + massRatio*massRatio));
  if (argument <= exp(beta*beta))
    dedx = 0.;
  else{
    dedx *= (log(argument)-beta*beta); // Bethe-Bloch [MeV/cm]
    dedx *= 1.E-3;  // in GeV/cm, hence 1.e-3
    assert(dedx>0);
  }

  double DE = step * dedx; //always positive
  double momLoss = sqrt(mom*mom+2.*sqrt(mom*mom+mass*mass)*DE+DE*DE) - mom; //always positive

  //in vacuum it can numerically happen that momLoss becomes a small negative number. A cut-off at 0.01 eV for momentum loss seems reasonable
  if(fabs(momLoss)<1.E-11)momLoss=1.E-11;
  return momLoss;
}


void GFMaterialEffects::noiseBetheBloch(const double& step,
                                        const double& mom,
                                        const int&    pdg,
                                        const double& matDensity,
                                        const double& matZ,
                                        const double& matA,
                                        const double& meanExcitationEnergy,
                                        TMatrixT<double>* noise){

  double charge, mass;
  getParticleParameters(pdg, charge, mass);

  const double beta = mom/sqrt(mass*mass+mom*mom);
  double dedx = 0.307075*matZ/matA*matDensity/(beta*beta)*charge*charge;

  //for numerical stability
  double gammaSquare = 1.-beta*beta;
  if(gammaSquare>1.E-10) gammaSquare = 1./gammaSquare;
  else gammaSquare = 1.E10;
  double gamma = sqrt(gammaSquare);

  double massRatio = me/mass;
  double argument = gammaSquare*beta*beta*me*1.E3*2./((1.E-6*meanExcitationEnergy) * sqrt(1+2*sqrt(gammaSquare)*massRatio + massRatio*massRatio));
  if (argument <= exp(beta*beta))
    dedx = 0.;
  else{
    dedx *= (log(argument)-beta*beta); // Bethe-Bloch [MeV/cm]
    dedx *= 1.E-3;  // in GeV/cm, hence 1.e-3
    assert(dedx>0);
  }

  // ENERGY LOSS FLUCTUATIONS; calculate sigma^2(E);
  double sigma2E = 0.;
  double zeta  = 153.4E3 * charge*charge/(beta*beta) * matZ/matA * matDensity * step; // eV
  double Emax  = 2.E9*me*beta*beta*gammaSquare / (1. + 2.*gamma*me/mass + (me/mass)*(me/mass) ); // eV
  double kappa = zeta/Emax;

  if (kappa > 0.01) { // Vavilov-Gaussian regime
    sigma2E += zeta*Emax*(1.-beta*beta/2.);  // eV^2
  }
  else { // Urban/Landau approximation
    double alpha = 0.996;
    double sigmaalpha = 15.76;
    // calculate number of collisions Nc
    double I = 16. * pow(matZ, 0.9); // eV
    double f2 = 0.;
    if (matZ > 2.) f2 = 2./matZ;
    double f1 = 1. - f2;
    double e2 = 10.*matZ*matZ; // eV
    double e1 = pow( (I/pow(e2,f2)), 1./f1);  // eV

    double mbbgg2 = 2.E9*mass*beta*beta*gammaSquare; // eV
    double Sigma1 = dedx*1.0E9 * f1/e1 * (log(mbbgg2 / e1) - beta*beta) / (log(mbbgg2 / I) - beta*beta) * 0.6; // 1/cm
    double Sigma2 = dedx*1.0E9 * f2/e2 * (log(mbbgg2 / e2) - beta*beta) / (log(mbbgg2 / I) - beta*beta) * 0.6; // 1/cm
    double Sigma3 = dedx*1.0E9 * Emax / ( I*(Emax+I)*log((Emax+I)/I) ) * 0.4; // 1/cm

    double Nc = (Sigma1 + Sigma2 + Sigma3)*step;

    if (Nc > 50.) { // truncated Landau distribution
      // calculate sigmaalpha  (see GEANT3 manual W5013)
      double RLAMED = -0.422784 - beta*beta - log(zeta/Emax);
      double RLAMAX =  0.60715 + 1.1934*RLAMED +(0.67794 + 0.052382*RLAMED)*exp(0.94753+0.74442*RLAMED);
      // from lambda max to sigmaalpha=sigma (empirical polynomial)
      if(RLAMAX <= 1010.) {
         sigmaalpha =  1.975560
                      +9.898841e-02 *RLAMAX
                      -2.828670e-04 *RLAMAX*RLAMAX
                      +5.345406e-07 *pow(RLAMAX,3.)
                      -4.942035e-10 *pow(RLAMAX,4.)
                      +1.729807e-13 *pow(RLAMAX,5.);
      }
      else { sigmaalpha = 1.871887E+01 + 1.296254E-02 *RLAMAX; }
      // alpha=54.6  corresponds to a 0.9996 maximum cut
      if(sigmaalpha > 54.6) sigmaalpha=54.6;
      sigma2E += sigmaalpha*sigmaalpha * zeta*zeta;  // eV^2
    }
    else { // Urban model
      double Ealpha  = I / (1.-(alpha*Emax/(Emax+I)));   // eV
      double meanE32 = I*(Emax+I)/Emax * (Ealpha-I);     // eV^2
      sigma2E += step * (Sigma1*e1*e1 + Sigma2*e2*e2 + Sigma3*meanE32); // eV^2
    }
  }

  sigma2E*=1.E-18; // eV -> GeV

  // update noise matrix
  (*noise)[6][6] += (mom*mom+mass*mass)/pow(mom,6.)*sigma2E;
}


void GFMaterialEffects::noiseCoulomb(const double& step,
                                     const double& mom,
                                     const int&    pdg,
                                     const double& matZ,
                                     const double& radiationLength,
                                           TMatrixT<double>* noise,
                                     const TMatrixT<double>* jacobian,
                                     const TVector3* directionBefore,
                                     const TVector3* directionAfter){

  double charge, mass;
  getParticleParameters(pdg, charge, mass);

  const double beta = mom/sqrt(mass*mass+mom*mom);

  // MULTIPLE SCATTERING; calculate sigma^2
  // PANDA report PV/01-07 eq(43); linear in step length
  double sigma2 = 225.E-6/(beta*beta*mom*mom) * step/radiationLength * matZ/(matZ+1) * log(159.*pow(matZ,-1./3.))/log(287.*pow(matZ,-0.5)); // sigma^2 = 225E-6/mom^2 * XX0/beta^2 * Z/(Z+1) * ln(159*Z^(-1/3))/ln(287*Z^(-1/2)

  // noiseBefore
    TMatrixT<double> noiseBefore(7,7);

    // calculate euler angles theta, psi (so that directionBefore' points in z' direction)
    double psi = 0;
    if (fabs((*directionBefore)[1]) < 1E-14) {  // numerical case: arctan(+-inf)=+-PI/2
      if ((*directionBefore)[0] >= 0.) psi = M_PI/2.;
      else psi = 3.*M_PI/2.;
    }
    else {
      if ((*directionBefore)[1] > 0.) psi = M_PI - atan((*directionBefore)[0]/(*directionBefore)[1]);
      else psi = -atan((*directionBefore)[0]/(*directionBefore)[1]);
    }
    // cache sin and cos
    double sintheta = sqrt(1-(*directionBefore)[2]*(*directionBefore)[2]);  // theta = arccos(directionBefore[2])
    double costheta = (*directionBefore)[2];
    double sinpsi = sin(psi);
    double cospsi = cos(psi);

    // calculate NoiseBefore Matrix R M R^T
    double noiseBefore34 =  sigma2 * cospsi * sinpsi * sintheta*sintheta; // noiseBefore_ij = noiseBefore_ji
    double noiseBefore35 = -sigma2 * costheta * sinpsi * sintheta;
    double noiseBefore45 =  sigma2 * costheta * cospsi * sintheta;

    noiseBefore[3][3] = sigma2 * (cospsi*cospsi + costheta*costheta - costheta*costheta * cospsi*cospsi);
    noiseBefore[4][3] = noiseBefore34;
    noiseBefore[5][3] = noiseBefore35;

    noiseBefore[3][4] = noiseBefore34;
    noiseBefore[4][4] = sigma2 * (sinpsi*sinpsi + costheta*costheta * cospsi*cospsi);
    noiseBefore[5][4] = noiseBefore45;

    noiseBefore[3][5] = noiseBefore35;
    noiseBefore[4][5] = noiseBefore45;
    noiseBefore[5][5] = sigma2 * sintheta*sintheta;

    TMatrixT<double> jacobianT(7,7);
    jacobianT = (*jacobian);
    jacobianT.T();

    noiseBefore = jacobianT*noiseBefore*(*jacobian); //propagate

  // noiseAfter
    TMatrixT<double> noiseAfter(7,7);

    // calculate euler angles theta, psi (so that A' points in z' direction)
    psi = 0;
    if (fabs((*directionAfter)[1]) < 1E-14) {  // numerical case: arctan(+-inf)=+-PI/2
      if ((*directionAfter)[0] >= 0.) psi = M_PI/2.;
      else psi = 3.*M_PI/2.;
    }
    else {
      if ((*directionAfter)[1] > 0.) psi = M_PI - atan((*directionAfter)[0]/(*directionAfter)[1]);
      else psi = -atan((*directionAfter)[0]/(*directionAfter)[1]);
    }
    // cache sin and cos
    sintheta = sqrt(1-(*directionAfter)[2]*(*directionAfter)[2]);  // theta = arccos(directionAfter[2])
    costheta = (*directionAfter)[2];
    sinpsi = sin(psi);
    cospsi = cos(psi);

    // calculate NoiseAfter Matrix R M R^T
    double noiseAfter34 =  sigma2 * cospsi * sinpsi * sintheta*sintheta; // noiseAfter_ij = noiseAfter_ji
    double noiseAfter35 = -sigma2 * costheta * sinpsi * sintheta;
    double noiseAfter45 =  sigma2 * costheta * cospsi * sintheta;

    noiseAfter[3][3] = sigma2 * (cospsi*cospsi + costheta*costheta - costheta*costheta * cospsi*cospsi);
    noiseAfter[4][3] = noiseAfter34;
    noiseAfter[5][3] = noiseAfter35;

    noiseAfter[3][4] = noiseAfter34;
    noiseAfter[4][4] = sigma2 * (sinpsi*sinpsi + costheta*costheta * cospsi*cospsi);
    noiseAfter[5][4] = noiseAfter45;

    noiseAfter[3][5] = noiseAfter35;
    noiseAfter[4][5] = noiseAfter45;
    noiseAfter[5][5] = sigma2 * sintheta*sintheta;

  //calculate mean of noiseBefore and noiseAfter and update noise
    (*noise) += 0.5*noiseBefore + 0.5*noiseAfter;

}


double GFMaterialEffects::energyLossBrems(const double& step,
                                     const double& mom,
                                     const int&    pdg,
                                     const double& matDensity,
                                     const double& matZ,
                                     const double& matA,
                                     const double& radiationLength){

  if (fabs(pdg)!=11) return 0; // only for electrons and positrons

  #if !defined(BETHE)
    static const double C[101]={ 0.0,-0.960613E-01, 0.631029E-01,-0.142819E-01, 0.150437E-02,-0.733286E-04, 0.131404E-05, 0.859343E-01,-0.529023E-01, 0.131899E-01,-0.159201E-02, 0.926958E-04,-0.208439E-05,-0.684096E+01, 0.370364E+01,-0.786752E+00, 0.822670E-01,-0.424710E-02, 0.867980E-04,-0.200856E+01, 0.129573E+01,-0.306533E+00, 0.343682E-01,-0.185931E-02, 0.392432E-04, 0.127538E+01,-0.515705E+00, 0.820644E-01,-0.641997E-02, 0.245913E-03,-0.365789E-05, 0.115792E+00,-0.463143E-01, 0.725442E-02,-0.556266E-03, 0.208049E-04,-0.300895E-06,-0.271082E-01, 0.173949E-01,-0.452531E-02, 0.569405E-03,-0.344856E-04, 0.803964E-06, 0.419855E-02,-0.277188E-02, 0.737658E-03,-0.939463E-04, 0.569748E-05,-0.131737E-06,-0.318752E-03, 0.215144E-03,-0.579787E-04, 0.737972E-05,-0.441485E-06, 0.994726E-08, 0.938233E-05,-0.651642E-05, 0.177303E-05,-0.224680E-06, 0.132080E-07,-0.288593E-09,-0.245667E-03, 0.833406E-04,-0.129217E-04, 0.915099E-06,-0.247179E-07, 0.147696E-03,-0.498793E-04, 0.402375E-05, 0.989281E-07,-0.133378E-07,-0.737702E-02, 0.333057E-02,-0.553141E-03, 0.402464E-04,-0.107977E-05,-0.641533E-02, 0.290113E-02,-0.477641E-03, 0.342008E-04,-0.900582E-06, 0.574303E-05, 0.908521E-04,-0.256900E-04, 0.239921E-05,-0.741271E-07,-0.341260E-04, 0.971711E-05,-0.172031E-06,-0.119455E-06, 0.704166E-08, 0.341740E-05,-0.775867E-06,-0.653231E-07, 0.225605E-07,-0.114860E-08,-0.119391E-06, 0.194885E-07, 0.588959E-08,-0.127589E-08, 0.608247E-10};
    static const double xi=2.51, beta=0.99, vl=0.00004;
  #endif
  #if defined(BETHE) // no MIGDAL corrections
    static const double C[101]={ 0.0, 0.834459E-02, 0.443979E-02,-0.101420E-02, 0.963240E-04,-0.409769E-05, 0.642589E-07, 0.464473E-02,-0.290378E-02, 0.547457E-03,-0.426949E-04, 0.137760E-05,-0.131050E-07,-0.547866E-02, 0.156218E-02,-0.167352E-03, 0.101026E-04,-0.427518E-06, 0.949555E-08,-0.406862E-02, 0.208317E-02,-0.374766E-03, 0.317610E-04,-0.130533E-05, 0.211051E-07, 0.158941E-02,-0.385362E-03, 0.315564E-04,-0.734968E-06,-0.230387E-07, 0.971174E-09, 0.467219E-03,-0.154047E-03, 0.202400E-04,-0.132438E-05, 0.431474E-07,-0.559750E-09,-0.220958E-02, 0.100698E-02,-0.596464E-04,-0.124653E-04, 0.142999E-05,-0.394378E-07, 0.477447E-03,-0.184952E-03,-0.152614E-04, 0.848418E-05,-0.736136E-06, 0.190192E-07,-0.552930E-04, 0.209858E-04, 0.290001E-05,-0.133254E-05, 0.116971E-06,-0.309716E-08, 0.212117E-05,-0.103884E-05,-0.110912E-06, 0.655143E-07,-0.613013E-08, 0.169207E-09, 0.301125E-04,-0.461920E-04, 0.871485E-05,-0.622331E-06, 0.151800E-07,-0.478023E-04, 0.247530E-04,-0.381763E-05, 0.232819E-06,-0.494487E-08,-0.336230E-04, 0.223822E-04,-0.384583E-05, 0.252867E-06,-0.572599E-08, 0.105335E-04,-0.567074E-06,-0.216564E-06, 0.237268E-07,-0.658131E-09, 0.282025E-05,-0.671965E-06, 0.565858E-07,-0.193843E-08, 0.211839E-10, 0.157544E-04,-0.304104E-05,-0.624410E-06, 0.120124E-06,-0.457445E-08,-0.188222E-05,-0.407118E-06, 0.375106E-06,-0.466881E-07, 0.158312E-08, 0.945037E-07, 0.564718E-07,-0.319231E-07, 0.371926E-08,-0.123111E-09};
    static const double xi=2.10, beta=1.00, vl=0.001;
  #endif

  double BCUT=10000.; // energy up to which soft bremsstrahlung energy loss is calculated

  double charge, mass;
  getParticleParameters(pdg, charge, mass);

  double THIGH=100., CHIGH=50.;
  double dedxBrems=0.;

  if(BCUT>0.){
    double T, kc;

    if(BCUT>=mom) BCUT=mom; // confine BCUT to mom

    // T=mom,  confined to THIGH
    // kc=BCUT, confined to CHIGH ??
    if(mom>=THIGH) {
      T=THIGH;
      if(BCUT>=THIGH) kc=CHIGH;
      else kc=BCUT;
    }
    else {
      T=mom;
      kc=BCUT;
    }

    double E=T+me; // total electron energy
    if(BCUT>T) kc=T;

    double X=log(T/me);
    double Y=log(kc/(E*vl));

    double XX;
    int    K;
    double S=0., YY=1.;

    for (unsigned int I=1; I<=2; ++I) {
      XX=1.;
      for (unsigned int J=1; J<=6; ++J) {
        K=6*I+J-6;
        S=S+C[K]*XX*YY;
        XX=XX*X;
      }
      YY=YY*Y;
    }

    for (unsigned int I=3; I<=6; ++I) {
       XX=1.;
       for (unsigned int J=1; J<=6; ++J) {
          K=6*I+J-6;
          if(Y<=0.) S=S+C[K]*XX*YY;
          else      S=S+C[K+24]*XX*YY;
          XX=XX*X;
       }
       YY=YY*Y;
    }

    double SS=0.;
    YY=1.;

    for (unsigned int I=1; I<=2; ++I) {
      XX=1.;
      for (unsigned int J=1; J<=5; ++J) {
        K=5*I+J+55;
        SS=SS+C[K]*XX*YY;
        XX=XX*X;
      }
      YY=YY*Y;
    }

    for (unsigned int I=3; I<=5; ++I) {
      XX=1.;
      for (unsigned int J=1; J<=5; ++J) {
        K=5*I+J+55;
        if(Y<=0.) SS=SS+C[K]*XX*YY;
        else      SS=SS+C[K+15]*XX*YY;
        XX=XX*X;
      }
      YY=YY*Y;
    }

    S=S+matZ*SS;

    if(S>0.){
      double CORR=1.;
      #if !defined(BETHE)
        CORR=1./(1.+0.805485E-10*matDensity*matZ*E*E/(matA*kc*kc)); // MIGDAL correction factor
      #endif

      double FAC=matZ*(matZ+xi)*E*E * pow((kc*CORR/T),beta) / (E+me);
      if(FAC<=0.) return 0.;
      dedxBrems=FAC*S;

      double RAT;

      if(mom>THIGH) {
        if(BCUT<THIGH) {
          RAT=BCUT/mom;
          S=(1.-0.5*RAT+2.*RAT*RAT/9.);
          RAT=BCUT/T;
          S=S/(1.-0.5*RAT+2.*RAT*RAT/9.);
        }
        else {
          RAT=BCUT/mom;
          S=BCUT*(1.-0.5*RAT+2.*RAT*RAT/9.);
          RAT=kc/T;
          S=S/(kc*(1.-0.5*RAT+2.*RAT*RAT/9.));
        }
        dedxBrems=dedxBrems*S; // GeV barn
      }

      dedxBrems = 0.60221367*matDensity*dedxBrems/matA; // energy loss dE/dx [GeV/cm]
    }
  }

  assert(dedxBrems>=0.);

  double factor=1.; // positron correction factor

  if (pdg==-11){
      static const double AA=7522100., A1=0.415, A3=0.0021, A5=0.00054;

      double ETA=0.;
      if(matZ>0.) {
        double X=log(AA*mom/matZ*matZ);
        if(X>-8.) {
          if(X>=+9.) ETA=1.;
          else {
             double W=A1*X+A3*pow(X,3.)+A5*pow(X,5.);
             ETA=0.5+atan(W)/M_PI;
          }
        }
      }

      double E0;

      if(ETA<0.0001) factor=1.E-10;
      else if (ETA>0.9999) factor=1.;
      else {
         E0=BCUT/mom;
         if(E0>1.) E0=1.;
         if(E0<1.E-8) factor=1.;
         else factor = ETA * ( 1.-pow(1.-E0, 1./ETA) ) / E0;
      }
  }

  double DE = step * factor*dedxBrems; //always positive
  double momLoss = sqrt(mom*mom+2.*sqrt(mom*mom+mass*mass)*DE+DE*DE) - mom; //always positive

  return momLoss;
}


void GFMaterialEffects::noiseBrems(const double& step,
                                   const double& mom,
                                   const int&    pdg,
                                   const double& radiationLength,
                                   TMatrixT<double>* noise){

  if (fabs(pdg)!=11) return; // only for electrons and positrons

  double charge, mass;
  getParticleParameters(pdg, charge, mass);

  double LX  = 1.442695*step/radiationLength;
  double S2B = mom*mom * ( 1./pow(3.,LX) - 1./pow(4.,LX) );
  double DEDXB  = pow(fabs(S2B),0.5);
  DEDXB = 1.2E9*DEDXB; //eV
  double sigma2E = DEDXB*DEDXB; //eV^2
  sigma2E*=1.E-18; // eV -> GeV

  (*noise)[6][6] += (mom*mom+mass*mass)/pow(mom,6.)*sigma2E;
}

//---------------------------------------------------------------------------------

ClassImp(GFMaterialEffects)


/*
Reference for elemental mean excitation energies at:
http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
*/

const int MeanExcEnergy_NELEMENTS = 93; // 0 = vacuum, 1 = hydrogen, 92 = uranium
const float MeanExcEnergy_vals[] = {1.e15, 19.2, 41.8, 40.0, 63.7, 76.0, 78., 82.0, 95.0, 115.0, 137.0, 149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0, 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0};

float MeanExcEnergy_get(int Z){
    assert(Z>=0&&Z<MeanExcEnergy_NELEMENTS);
    return MeanExcEnergy_vals[Z];
}

float MeanExcEnergy_get(TGeoMaterial* mat){
  if(mat->IsMixture()){
    double logMEE = 0.;
    double denom  = 0.;
    TGeoMixture *mix = (TGeoMixture*)mat;
    for(int i=0;i<mix->GetNelements();++i){
      int index = int(floor((mix->GetZmixt())[i]));
      logMEE += 1./(mix->GetAmixt())[i]*(mix->GetWmixt())[i]*(mix->GetZmixt())[i]*log(MeanExcEnergy_get(index));
      denom  += (mix->GetWmixt())[i]*(mix->GetZmixt())[i]*1./(mix->GetAmixt())[i];
    }
    logMEE/=denom;
    return exp(logMEE);
  }
  else{ // not a mixture
    int index = int(floor(mat->GetZ()));
    return MeanExcEnergy_get(index);
  }
}
