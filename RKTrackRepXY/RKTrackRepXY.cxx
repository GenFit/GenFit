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
#include"RKTrackRepXY.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include"assert.h"
#include"math.h"
#include"TMath.h"
#include"TGeoManager.h"
#include"TDatabasePDG.h"
#include"MeanExcEnergy.h"
#include"energyLoss.h"
#include"GFException.h"
#include"GFFieldManager.h"

RKTrackRepXY::~RKTrackRepXY(){
}

RKTrackRepXY::RKTrackRepXY() : GFAbsTrackRep(5),fPdg(0),fMass(0.),fDirection(true){
}

RKTrackRepXY::RKTrackRepXY(const TVector3& pos,
		       const TVector3& mom,
		       const TVector3& poserr,
		       const TVector3& momerr,
		       //const double& q,
		       const int& PDGCode) :
  GFAbsTrackRep(5),fDirection(true){
  //sets charge q also
  setPDG(PDGCode);
  TVector3 orig(0.,0.,pos.Z());
  static const TVector3 u(1.,0.,0.);
  static const TVector3 v(0.,1.,0.);
  fRefPlane.setO(orig);
  fRefPlane.setU(u);
  fRefPlane.setV(v);

  fState[0][0]=pos.X();
  fState[1][0]=pos.Y();
  fState[2][0]=mom.X()/mom.Z();
  fState[3][0]=mom.Y()/mom.Z();
  fState[4][0]=fCharge/mom.Mag();
  
  fCov[0][0] = poserr.X()*poserr.X();
  fCov[1][1] = poserr.Y()*poserr.Y();
  fCov[2][2] = 1./(mom.Z()*mom.Z())*
    (momerr.X()*momerr.X()+
     mom.X()*mom.X()/(mom.Z()*mom.Z())*momerr.Z()*momerr.Z());
  fCov[3][3] = 1./(mom.Z()*mom.Z())*
    (momerr.Y()*momerr.Y()+
     mom.Y()*mom.Y()/(mom.Z()*mom.Z())*momerr.Z()*momerr.Z());
  assert(fabs(fCharge)>1.E-3);//charged particle
  fCov[4][4] = fCharge*fCharge/pow(mom.Mag(),6.)*
    (mom.X()*mom.X()*momerr.X()*momerr.X()+
     mom.Y()*mom.Y()*momerr.Y()*momerr.Y()+
     mom.Z()*mom.Z()*momerr.Z()*momerr.Z());
}

void RKTrackRepXY::setPDG(int i){
  fPdg = i;
  TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(fPdg);
  if(part == 0){
    std::cerr << "RKTrackRepXY::setPDG particle " << i 
	      << " not known to TDatabasePDG -> abort" << std::endl;
    exit(1);
  }
  fMass = part->Mass();
  fCharge = part->Charge()/3.;
  //  std::cout << part->GetName() << " " << charge << std::endl;
}

TVector3 RKTrackRepXY::getPos(const GFDetPlane& pl){
  if(pl!=fRefPlane){
    TMatrixT<double> s(5,1);
    TMatrixT<double> c(5,5);
    extrapolate(pl,s,c);
    return pl.getO()+s[0][0]*pl.getU()+s[1][0]*pl.getV();
  }
  return fRefPlane.getO()+fState[0][0]*fRefPlane.getU()+fState[1][0]*fRefPlane.getV();
}
TVector3 RKTrackRepXY::getMom(const GFDetPlane& pl){
  if(pl!=fRefPlane){
    TMatrixT<double> s(5,1);
    TMatrixT<double> c(5,5);
    extrapolate(pl,s,c);
    TVector3 retVal(0.,0.,1.);
    retVal += s[2][0]*pl.getU()+s[3][0]*pl.getV();
    retVal.SetMag(fabs(1./s[4][0]));
    return retVal;
  }
  else{
    TVector3 retVal(0.,0.,1.);
    retVal += fState[2][0]*fRefPlane.getU()+fState[3][0]*fRefPlane.getV();
    retVal.SetMag(fabs(1./fState[4][0]));
    return retVal;
  }
}
void RKTrackRepXY::getPosMom(const GFDetPlane& pl,TVector3& pos,
			   TVector3& mom){
  mom.SetXYZ(fState[2][0],fState[3][0],1.);
  mom.SetMag(fabs(1./fState[4][0]));
  pos.SetXYZ(fState[0][0],fState[1][0],fRefPlane.getO().Z());
  if(pl!=fRefPlane){
    TMatrixT<double> s(5,1);
    TMatrixT<double> c(5,5);
    extrapolate(pl,s,c);
    mom.SetXYZ(s[2][0],s[3][0],1.);
    mom.SetMag(fabs(1./s[4][0]));
    pos.SetXYZ(s[0][0],s[1][0],pl.getO().Z());
  }
  else{
    mom.SetXYZ(fState[2][0],fState[3][0],1.);
    mom.SetMag(fabs(1./fState[4][0]));
    pos.SetXYZ(fState[0][0],fState[1][0],fRefPlane.getO().Z());
  }
}

double RKTrackRepXY::getStep(const double& zFinal, const double& distance, const TVector3& pos, const TVector3& mom,double& XX0, double& dP){
  double dist = distance;
  double dz = zFinal-pos.Z();
  TVector3 dir = mom;
  dir.SetMag(1.);
  if(dir.Z() > 0){//particle flying downstream
    if(dz<0.){//backward extrap
      dir *= -1.;
    }
  }
  else{//flying upstream
    if(dz>0.){//backward extrap
      dir *= -1.;
    }
  }
  

  gGeoManager->InitTrack(pos.X(),pos.Y(),pos.Z(),dir.X(),dir.Y(),dir.Z());

  XX0 = 0.;  
  double X(0.);
  dP = 0.;
  double p = mom.Mag();
  while(dist>1.e-3){//10 micron
    assert(gGeoManager->GetCurrentVolume()->GetMedium()!=NULL);
    TGeoMaterial * mat = gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial();
    double radLen = mat->GetRadLen();
    //mat->Print();
    //dont calculate dedx for Z==0, i.e. vacuum
    double dedx = 0.;
    if(mat->GetZ()>1.E-3){
      dedx = energyLoss(p/sqrt(fMass*fMass+p*p),//beta
			getCharge(),//particle charge
			mat->GetDensity(),//material density
			mat->GetZ()/mat->GetA(),//material Z/A
			MeanExcEnergy::get(mat)//mean exc energy
			);
    }
    gGeoManager->FindNextBoundaryAndStep(dist);
    double step = gGeoManager->GetStep();
    double DE = step * dedx;
    double DP = p - sqrt(p*p-2*sqrt(p*p+fMass*fMass)*DE+DE*DE);
    static const double maxPloss = .01;
    if(dP + DP > p*maxPloss){
      double fraction = (p*maxPloss-dP)/DP;
      assert(fraction>0.&&fraction<1.);
      dP+=fraction*DP;
      X+=fraction*step;
      XX0+=fraction*step/radLen;
      break;
    }
    dP += DP;
    XX0 += step/radLen;
    X += step;
    dist -= step;
  }
  return (pos+X*dir).Z();
}

void RKTrackRepXY::extrapolateToPoint(const TVector3& pos,
				   TVector3& poca,
				   TVector3& dirInPoca){
  GFDetPlane d;
  d.setU(1.,0.,0.);
  d.setV(0.,1.,0.);
  d.setO(0.,0.,pos.Z());
  TMatrixT<double> s(5,1);
  TMatrixT<double> c(5,5);
  extrapolate(d,s,c);
  poca.SetXYZ(s[0][0],s[1][0],pos.Z());
  dirInPoca.SetXYZ(0.,0.,1.);
}
  
void RKTrackRepXY::extrapolateToLine(const TVector3& point1,
				   const TVector3& point2,
				   TVector3& poca,
				   TVector3& dirInPoca,
				   TVector3& poca_onwire){
  TVector3 line = point1-point2;
  line.SetMag(1.);
  if(asin(line.Z())/TMath::Pi()*180.>1.){
    GFException exc("asin(line.Z())/TMath::Pi()*180.>1. - wire needs to be perpendicular to z-axis",__LINE__,__FILE__);
    throw exc;
  }
  double z = (0.5*(point1+point2)).Z();
  GFDetPlane d;
  d.setU(1.,0.,0.);
  d.setV(0.,1.,0.);
  d.setO(0.,0.,z);
  TMatrixT<double> s(5,1);
  TMatrixT<double> c(5,5);
  extrapolate(d,s,c);
  poca.SetXYZ(s[0][0],s[1][0],z);
  dirInPoca.SetXYZ(0.,0.,1.);
  double t = 1./((point2-point1)*(point2-point1))*(poca*(point2-point1)+point1*point1-point1*point2);
  poca_onwire = point1+t*(point2-point1);
}


double RKTrackRepXY::extrapolate(const GFDetPlane& pl, 
			       TMatrixT<double>& statePred,
			       TMatrixT<double>& covPred){
  //perpendicularity cut for plane being perp to the z-axis
  static const double cosPerpCut = cos(1/180.*TMath::Pi());
  if(fabs(pl.getNormal().Z())<cosPerpCut){
    GFException exc("RKTrackRepXY only works for GFDetPlane perpendicular to beam axis - fabs(pl.getNormal().Z())<cosPerpCut",__LINE__,__FILE__);
    //exc.setFatal();
    throw exc;
  }

  TVector3 O = fRefPlane.getO();
  TVector3 U = fRefPlane.getU();
  TVector3 V = fRefPlane.getV();

  TMatrixT<double> rot(5,5);
  double th;
  if(U.Y()>=0.) th = acos(U.X());
  else th = TMath::Pi() + acos(-1.*U.X());

  rot[0][0] = cos(th);
  rot[1][1] = cos(th);
  rot[2][2] = cos(th);
  rot[3][3] = cos(th);
  rot[0][1] = -1.*sin(th);
  rot[1][0] = sin(th);
  rot[2][3] = -1.*sin(th);
  rot[3][2] = sin(th);
  rot[4][4] = 1.;

  TMatrixT<double> trans(5,1);
  trans[0][0] = O.X();
  trans[1][0] = O.Y();

  TMatrixT<double> stateExtrapIn = rot*fState;
  stateExtrapIn += trans;

  TMatrixT<double> rotT(rot);
  rotT.T();
  TMatrixT<double> rotCov = rot * ( fCov * rotT );

  //==============================
  TMatrixT<double> covExtrapIn = cov25to15(rotCov);

  TMatrixT<double> covExtrapOut(15,1);
  TMatrixT<double> stateExtrapOut(5,1);
  double zFrom(fRefPlane.getO().Z());
  double zExtrap(0.);
  double distance(0.);
  int counter(0);
  while(true){
    double dz = pl.getO().Z()-zFrom;
    double dist = sqrt(dz*dz*(1.+stateExtrapIn[2][0]*stateExtrapIn[2][0]+stateExtrapIn[3][0]*stateExtrapIn[3][0]));
    
    TVector3 pos(stateExtrapIn[0][0],stateExtrapIn[1][0],zFrom);
    TVector3 mom(stateExtrapIn[2][0],stateExtrapIn[3][0],1.);
    mom.SetMag(fabs(1./stateExtrapIn[4][0]));
    
    double XX0,dP;
    zExtrap = getStep(pl.getO().Z(), dist,  pos,  mom, XX0,  dP);
    
    bool breakFlag(false);
    static const double CUT = 1.e-3;
    if(fabs(zExtrap-pl.getO().Z())<CUT){
      zExtrap = pl.getO().Z();
      breakFlag = true;
    }

    if(dz>0.) {
      assert(zExtrap<pl.getO().Z()+CUT);
    }
    else{
      assert(zExtrap>pl.getO().Z()-CUT);
    }

    double zFinal(-1.E300);
    double thisDistance = this->Extrap(zExtrap,zFrom,zFinal,stateExtrapIn,stateExtrapOut,covExtrapIn,covExtrapOut);
    //correct radiation length for helix path length as compared to straight line
    distance+=thisDistance;
    XX0*=thisDistance/dist;
    //add multiple scattering to covariance matrix
    addNoise(XX0,stateExtrapIn,covExtrapOut);
    //  std::cout << "dP/P " << dP/mom.Mag() << std::endl;
    if(fDirection)  stateExtrapOut[4][0] = getCharge()/(mom.Mag()-dP);
    else stateExtrapOut[4][0] = getCharge()/(mom.Mag()+dP);
    stateExtrapIn = stateExtrapOut;
    covExtrapIn = covExtrapOut;
    zFrom = zExtrap;
    if(breakFlag) break;
  }
  //printf("in x,y,z: %f,%f,%f\n",stateExtrapOut[0][0],
  // stateExtrapOut[1][0],pl.getO().Z());
  //printf("mom in x,y,z:\n");
  //TVector3 temp(stateExtrapOut[2][0],stateExtrapOut[3][0],1.);
  //temp.SetMag(fabs(1./stateExtrapOut[4][0]));
  //temp.Print();

  //target plane
  //all variables already defined earlier at beginning
  O = pl.getO();
  U = pl.getU();
  V = pl.getV();

  if(U.Y()>=0.) th = acos(U.X());
  else th = TMath::Pi() + acos(-1.*U.X());
  //  std::cout << "THETA " << th/TMath::Pi() * 180. << std::endl;
  rot[0][0] = cos(th);
  rot[1][1] = cos(th);
  rot[2][2] = cos(th);
  rot[3][3] = cos(th);
  rot[0][1] = -1.*sin(th);
  rot[1][0] = sin(th);
  rot[2][3] = -1.*sin(th);
  rot[3][2] = sin(th);
  rot[4][4] = 1.;

  trans[0][0] = O.X();
  trans[1][0] = O.Y();

  rotT = rot;
  rotT.T();

  stateExtrapOut -= trans;
  statePred = rotT*stateExtrapOut;

  rotCov = cov15to25(covExtrapOut);
  covPred = rotT*rotCov*rot;;

  //pl.Print();
  //statePred.Print();
  //cov15to25(covExtrapOut).Print();
  //covPred.Print();  
  return distance;

}

/*
  virtual void extrapolateToPoint(const TVector3& point,
				 TVector3& poca,
				 TVector3& dirInPoca); 

  virtual void extrapolateToLine(const TVector3& point1, 
				 const TVector3& point2,
				 TVector3& poca,
				 TVector3& dirInPoca,
				 TVector3& poca_onwire);

 */
void RKTrackRepXY::addNoise(double len,const TMatrixT<double>& s, TMatrixT<double>& cov15){
  if(s[4][0] == 0.) return; // momentum not known. Do nothing.
  
  // Lynch and Dahl aproximation for Sigma(Theta_proj) of mult. scatt.
  double SigTheta = 0.0136*fabs(s[4][0]) * sqrt(len) * (1.+0.038*log(len));
 
  // Noise matrix calculation (NIM A329 (1993) 493-500)
  // Transverse displacement of the track is ignored.
  double p3 = s[2][0];
  double p4 = s[3][0];

  double p3p3 = SigTheta*SigTheta * (1 + p3*p3) * (1 + p3*p3 + p4*p4);
  double p4p4 = SigTheta*SigTheta * (1 + p4*p4) * (1 + p3*p3 + p4*p4);
  double p3p4 = SigTheta*SigTheta * p3*p4       * (1 + p3*p3 + p4*p4);

  cov15[5][0] = cov15[5][0]+p3p3;
  cov15[8][0] = cov15[8][0]+p3p4;
  cov15[9][0] = cov15[9][0]+p4p4;
}


TMatrixT<double> RKTrackRepXY::cov15to25(const TMatrixT<double>& cov15) const{
  TMatrixT<double> retVal(5,5);
  assert(cov15.GetNrows()==15);
  for(int j=0;j<5;++j){
    for(int i=j;i<5;++i){
      retVal[i][j] = cov15[j+((i+1)*i)/2][0];
      if(i!=j) retVal[j][i] = retVal[i][j];
    }
  }
  return retVal;
}

TMatrixT<double> RKTrackRepXY::cov25to15(const TMatrixT<double>& cov25) const{
  TMatrixT<double> retVal(15,1);
  assert(cov25.GetNrows()==5 && cov25.GetNcols()==5);
  for(int j=0;j<5;++j){
    for(int i=j;i<5;++i){
      retVal[j+((i+1)*i)/2][0] = cov25[i][j];
    }
  }  
  return retVal;
}


bool RKTrackRepXY::RKutta (double* SU,double* VO, double& Path) const {


  const double EC     = .00014989626;                // 
  const double DLT    = .0002;                       // 
  const double DLT32  = DLT/32.;                     //
  const double Sstop  = .001;                        // Min. step 
  const double P3     = .33333333;                   //
  const double Smax   = 100.;                        // max. step allowed>0 
  const double Wmax   = 5500.;                       // max. way allowed
  const int    ND     = 42, ND1=ND-7;                //
  double*  R          = &VO[0];                      // Start coordinates
  double*  A          = &VO[3];                      // Start directions
  double SA[3]        = {0.,0.,0.};                  // Start directions derivatives 
  double  Pinv        = VO[6]*EC;                    // Invert momentum/2.
  double  Way         = 0.;                          // Total way of the trajectory
  int     Error       = 0;                           // Error of propogation
  //
  // Step estimation until surface
  //
  double Step,An,Dist,Dis,S,Sl=0;
  Dist=SU[3]-R[0]*SU[0]-R[1]*SU[1]-R[2]*SU[2]; An=A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];
  if (An==0) return(false); Step=Dist/An; 
  if(fabs(Step)>Wmax) {
    std::cout<<"PaAlgo::RKutta ==> Too long extrapolation requested : "<<Step<<" cm !"<<std::endl;
    std::cout<<"X = "<<R[0]<<" Y = "<<R[1]<<" Z = "<<R[2]
	<<"  COSx = "<<A[0]<<"  COSy = "<<A[1]<<"  COSz = "<<A[2]<<std::endl;
    std::cout<<"Destination  X = "<<SU[0]*SU[3]<<std::endl;
    return(false);
  }
  Step>Smax ? S=Smax : Step<-Smax ? S=-Smax : S=Step;
  //
  // Main cycle of Runge-Kutta method
  //
  //
  while(fabs(Step)>Sstop) {
    double H0[12],H1[12],H2[12],r[3];
    double S3=P3*S, S4=.25*S, PS2=Pinv*S; 
    //
    // First point
    //   

    r[0]=R[0]      ; r[1]=R[1]      ; r[2]=R[2]      ;  
    TVector3 pos(r[1],r[2],r[0]);
    TVector3 H0vect = GFFieldManager::getFieldVal(pos);
    H0[0]=PS2*H0vect.Z(); H0[1]=PS2*H0vect.X(); H0[2]=PS2*H0vect.Y(); 
    double A0=A[1]*H0[2]-A[2]*H0[1], B0=A[2]*H0[0]-A[0]*H0[2], C0=A[0]*H0[1]-A[1]*H0[0];
    double A2=A[0]+A0              , B2=A[1]+B0              , C2=A[2]+C0              ;
    double A1=A2+A[0]              , B1=B2+A[1]              , C1=C2+A[2]              ;
    //
    // Second point
    //
    r[0]+=A1*S4    ; r[1]+=B1*S4    ; r[2]+=C1*S4    ;   //setup.Field(r,H1);
    pos.SetXYZ(r[1],r[2],r[0]);
    TVector3 H1vect = GFFieldManager::getFieldVal(pos);
    H1[0]=H1vect.Z()*PS2; H1[1]=H1vect.X()*PS2;H1[2]=H1vect.Y()*PS2;
    double A3,B3,C3,A4,B4,C4,A5,B5,C5;
    A3 = B2*H1[2]-C2*H1[1]+A[0]; B3=C2*H1[0]-A2*H1[2]+A[1]; C3=A2*H1[1]-B2*H1[0]+A[2];
    A4 = B3*H1[2]-C3*H1[1]+A[0]; B4=C3*H1[0]-A3*H1[2]+A[1]; C4=A3*H1[1]-B3*H1[0]+A[2];
    A5 = A4-A[0]+A4            ; B5=B4-A[1]+B4            ; C5=C4-A[2]+C4            ;
    //
    // Last point
    //
    r[0]=R[0]+S*A4 ; r[1]=R[1]+S*B4 ; r[2]=R[2]+S*C4 ;  //setup.Field(r,H2);
    pos.SetXYZ(r[1],r[2],r[0]);
    TVector3 H2vect = GFFieldManager::getFieldVal(pos);
    H2[0]=H2vect.Z()*PS2; H2[1]=H2vect.X()*PS2;H2[2]=H2vect.Y()*PS2;
    double A6=B5*H2[2]-C5*H2[1], B6=C5*H2[0]-A5*H2[2], C6=A5*H2[1]-B5*H2[0];
    //
    // Test approximation quality on give step and possible step reduction
    //
    double EST;
    if((EST=fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4)))>DLT) {S*=.5; continue;}
    //
    // Derivatives of track parameters in last point
    //
    for(int i=7; i!=ND; i+=7) {
      double* dR = &VO[i];  double* dA = &VO[i+3]; 
      double dA0   = H0[ 2]*dA[1]-H0[ 1]*dA[2];              // dA0/dp
      double dB0   = H0[ 0]*dA[2]-H0[ 2]*dA[0];              // dB0/dp
      double dC0   = H0[ 1]*dA[0]-H0[ 0]*dA[1];              // dC0/dp
      if(i==ND1) {dA0+=A0; dB0+=B0; dC0+=C0;}
      double dA2   = dA0+dA[0]; 
      double dB2   = dB0+dA[1]; 
      double dC2   = dC0+dA[2];
      double dA3   = dA[0]+dB2*H1[2]-dC2*H1[1];              // dA3/dp
      double dB3   = dA[1]+dC2*H1[0]-dA2*H1[2];              // dB3/dp
      double dC3   = dA[2]+dA2*H1[1]-dB2*H1[0];              // dC3/dp
      if(i==ND1) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];}
      double dA4   = dA[0]+dB3*H1[2]-dC3*H1[1];              // dA4/dp
      double dB4   = dA[1]+dC3*H1[0]-dA3*H1[2];              // dB4/dp
      double dC4   = dA[2]+dA3*H1[1]-dB3*H1[0];              // dC4/dp
      if(i==ND1) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];}
      double dA5   = dA4+dA4-dA[0];
      double dB5   = dB4+dB4-dA[1];
      double dC5   = dC4+dC4-dA[2];
      double dA6   = dB5*H2[2]-dC5*H2[1];                    // dA6/dp
      double dB6   = dC5*H2[0]-dA5*H2[2];                    // dB6/dp
      double dC6   = dA5*H2[1]-dB5*H2[0];                    // dC6/dp
      if(i==ND1) {dA6+=A6; dB6+=B6; dC6+=C6;}
      dR[0]+=(dA2+dA3+dA4)*S3; dA[0] = (dA0+dA3+dA3+dA5+dA6)*P3;      
      dR[1]+=(dB2+dB3+dB4)*S3; dA[1] = (dB0+dB3+dB3+dB5+dB6)*P3; 
      dR[2]+=(dC2+dC3+dC4)*S3; dA[2] = (dC0+dC3+dC3+dC5+dC6)*P3;
    }
    if((Way+=fabs(S))>Wmax){ 
      std::cout<<"PaAlgo::RKutta ==> Trajectory is longer then length limit : "<<Way<<" cm !"
	  <<" P = "<<1./VO[6]<< " GeV"<<std::endl;
      return(false);
    }
    //
    // Track parameters in last point
    //   
    R[0]+=(A2+A3+A4)*S3; A[0]+=(SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]); 
    R[1]+=(B2+B3+B4)*S3; A[1]+=(SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);
    R[2]+=(C2+C3+C4)*S3; A[2]+=(SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]); Sl=S;
    double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    A[0]*=CBA; A[1]*=CBA; A[2]*=CBA;
    //
    // Step estimation until surface and test conditions for stop propogation
    //
    Dis = SU[3]-R[0]*SU[0]-R[1]*SU[1]-R[2]*SU[2]; An=A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];
    if (An==0 || (Dis*Dist>0 && fabs(Dis)>fabs(Dist))) {Error=1; Step=0; break;}   
    Step = Dis/An; Dist=Dis;
    //
    // Possible current step reduction
    //

    if(EST<DLT32                      ) S*=2.;
    if(S*Step<0. || fabs(S)>fabs(Step)) S=Step;

  } //end of main loop

  //
  // Output information preparation for main track parameteres
  //
  if(Sl!=0) Sl=1./Sl;
  A [0]+=(SA[0]*=Sl)*Step; A [1]+=(SA[1]*=Sl)*Step; A [2]+=(SA[2]*=Sl)*Step;
  double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
  VO[0]      = R[0]+Step*(A[0]-.5*Step*SA[0]);                            
  VO[1]      = R[1]+Step*(A[1]-.5*Step*SA[1]);
  VO[2]      = R[2]+Step*(A[2]-.5*Step*SA[2]);
  VO[3]      = A[0]*CBA;
  VO[4]      = A[1]*CBA;
  VO[5]      = A[2]*CBA;
  //
  // Output derivatives of track parameters preparation 
  //
  An = A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2]; An!=0 ? An=1./An : An = 0;
  for(int i=7; i!=ND; i+=7) {
    double* dR = &VO[i];  double* dA = &VO[i+3];
    S = (dR[0]*SU[0]+dR[1]*SU[1]+dR[2]*SU[2])*An;
    dR[0]-=S*A [0];  dR[1]-=S*A [1]; dR[2]-=S*A [2]; 
    dA[0]-=S*SA[0];  dA[1]-=S*SA[1]; dA[2]-=S*SA[2];
  }
  if(Error == 1){
    std::cout<<"PaAlgo::RKutta ==> Do not getting closer. Path = "
	<<Way<<" cm"<<"  P = "<<1./VO[6]<<" GeV"<<std::endl;
    return(false);
  }

  Path=fabs(Way);
  return(true);
}


double RKTrackRepXY::Extrap( double Z, const double& zFrom,double& zOut, const TMatrixT<double>& stateIn, TMatrixT<double>& stateOut, const TMatrixT<double>& cov15,TMatrixT<double>& covOut) const{

  // for historical reasons here GEANT3 coodinate system is used 
  // (X along the beam, Z - upward)

  stateOut.ResizeTo(5,1);
  covOut.ResizeTo(15,1);
  assert(cov15.GetNrows()==15);  

  const double RKuttaMinStep = 0.0100;  // Ringe Kutta minimal step (100 mic)
  double xlast = Z;

  double xfirst = zFrom;
  double s(0);

  if(fabs(xlast-xfirst) < RKuttaMinStep) { // very short step
    zOut = xlast;
    stateOut = stateIn;
    covOut = cov15;
    return fabs(xlast-xfirst);
  }

  //
  // Undefined momentum. Do straight line extrapolation.
  //
  if(stateIn[4][0] == 0.){
    double dx=xlast-xfirst;
    zOut = xlast;
    stateOut[0][0]=stateIn[0][0]+stateIn[2][0]*dx;
    stateOut[1][0]=stateIn[1][0]+stateIn[3][0]*dx;
    stateOut[2][0]=stateIn[2][0];
    stateOut[3][0]=stateIn[3][0];
    stateOut[4][0]=stateIn[4][0];
    
    // Cov matrix propagation
    covOut[0][0] = cov15[0][0] + cov15[3][0]*dx + (cov15[3][0] + cov15[5][0]*dx)*dx;
    covOut[1][0] = cov15[1][0] + cov15[6][0]*dx + (cov15[4][0] + cov15[8][0]*dx)*dx;
    covOut[2][0] = cov15[2][0] + cov15[7][0]*dx + (cov15[7][0] + cov15[9][0]*dx)*dx;
    covOut[3][0] = cov15[3][0] + cov15[5][0]*dx;
    covOut[4][0] = cov15[4][0] + cov15[8][0]*dx;
    covOut[5][0] = cov15[5][0];
    covOut[6][0] = cov15[6][0] + cov15[8][0]*dx;
    covOut[7][0] = cov15[7][0] + cov15[9][0]*dx;
    covOut[8][0] = cov15[8][0];
    covOut[9][0] = cov15[9][0];
    covOut[10][0]= cov15[10][0]+ cov15[12][0]*dx;
    covOut[11][0]= cov15[11][0]+ cov15[13][0]*dx;
    covOut[12][0]= cov15[12][0];
    covOut[13][0]= cov15[13][0];
    covOut[14][0]= cov15[14][0];
    
    return sqrt(dx*dx+
		(stateOut[0][0]-stateIn[0][0])*(stateOut[0][0]-stateIn[0][0])+
		(stateOut[1][0]-stateIn[1][0])*(stateOut[1][0]-stateIn[1][0])
		);
  }
    
  //
  // Runge-Kutta extrapolation
  //

  // cut on  P < 100 MeV
  if(fabs(1./stateIn[4][0]) < 0.100) {
    std::ostringstream ostr;
    ostr<<"RKTrackRepXY::Extrap() ==> Too low momentum for Runge-Kutta propagation: "
	<<fabs(1./stateIn[4][0]) <<" GeV";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  

  // Prepare surface description array
  double su[4];
  if(xlast < 0 ) su[0] = -1.;
  else           su[0] =  1.;
  su[1]=su[2]=0.;
  su[3] = fabs(xlast);

  // Prepare track parameters

  double x  = zFrom;
  double y  = stateIn[0][0];
  double z  = stateIn[1][0];
  double ax = 1./sqrt(stateIn[2][0]*stateIn[2][0] + stateIn[3][0]*stateIn[3][0] + 1);
  double ay = stateIn[2][0]*ax;
  double az = stateIn[3][0]*ax;
  double qP = stateIn[4][0];
  
  // Prepare Jacobian dP/dH where 
  //
  // P is track parameters vector (size 7), needed for RKutta routine,
  // in the direction cosines (ax,ay,az) representation: (x,y,x,ax,ay,az,q/P)  
  //
  // H is helix.
  //

  double P[42] = {
    x,      y,      z,         ax,          ay,         az,     qP,   // Parameters (P) 
    0,      1,      0,          0,           0,           0,     0,   // dP/dH[1]			  
    0,      0,      1,          0,           0,           0,     0,   // dP/dH[2]
    0,      0,      0,  -ay*ax*ax, ax-ay*ay*ax,   -ay*az*ax,     0,   // dP/dH[3]
    0,      0,      0,  -az*ax*ax,   -ay*az*ax, ax-az*az*ax,     0,   // dP/dH[4]
    0,      0,      0,          0,           0,           0,     qP   // dP/dH[5]*H[5]
  };

  //
  // ---- Do the propagation
  //
  if( ! this->RKutta(su,P,s) ) {
    GFException exc("RKTrackRepXY Runge Kutta propagation failed",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  //
  // ----
  //

  if(P[3] <= 0.) {
    std::ostringstream ostr;
    ostr<<"RKTrackRepXY::Extrap ==> Trajectory had turned back during extrapolation. P = "
	<<fabs(1./P[6])<<" GeV";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // Output helix parameters

  zOut = P[0]; // x
  stateOut[0][0] = P[1]; // y
  stateOut[1][0] = P[2]; // z
  stateOut[2][0] = P[4]/P[3];    // yp = ay/ax
  stateOut[3][0] = P[5]/P[3];    // zp = az/ax
  stateOut[4][0] = P[6];


  // Calculate helix transformation Jacobian (F = dHout/dHin)  

   
  double p3Ax = -P[4]/(P[3]*P[3]);
  double p3Ay = 1./P[3];
  double p3Az = 0.;
  double p4Ax = -P[5]/(P[3]*P[3]);
  double p4Ay = 0.;
  double p4Az = p3Ay;
  
  
  double f[5][5];

  f[0][0]= P[ 8];    f[1][0]= P[ 9];
  f[0][1]= P[15];    f[1][1]= P[16];
  f[0][2]= P[22];    f[1][2]= P[23];
  f[0][3]= P[29];    f[1][3]= P[30];
  f[0][4]= P[36]/qP; f[1][4]= P[37]/qP;
  
  f[2][0]= p3Ax*P[10]+p3Ay*P[11]+p3Az*P[12];     f[3][0]= p4Ax*P[10]+p4Ay*P[11]+p4Az*P[12];   
  f[2][1]= p3Ax*P[17]+p3Ay*P[18]+p3Az*P[19];     f[3][1]= p4Ax*P[17]+p4Ay*P[18]+p4Az*P[19];
  f[2][2]= p3Ax*P[24]+p3Ay*P[25]+p3Az*P[26];     f[3][2]= p4Ax*P[24]+p4Ay*P[25]+p4Az*P[26];
  f[2][3]= p3Ax*P[31]+p3Ay*P[32]+p3Az*P[33];     f[3][3]= p4Ax*P[31]+p4Ay*P[32]+p4Az*P[33];
  f[2][4]=(p3Ax*P[38]+p3Ay*P[39]+p3Az*P[40])/qP; f[3][4]=(p4Ax*P[38]+p4Ay*P[39]+p4Az*P[40])/qP;
  
  f[4][0]=f[4][1]=f[4][2]=f[4][3]=0; f[4][4]=1;


  // Propagate Cov matrix (F*Cov*F.t)

  double w[5];

  w[0]=cov15[0][0] * f[0][0] + cov15[1][0] * f[0][1] + cov15[3][0] * f[0][2] + cov15[6][0] * f[0][3] + cov15[10][0]*f[0][4];
  w[1]=cov15[1][0] * f[0][0] + cov15[2][0] * f[0][1] + cov15[4][0] * f[0][2] + cov15[7][0] * f[0][3] + cov15[11][0]*f[0][4];
  w[2]=cov15[3][0] * f[0][0] + cov15[4][0] * f[0][1] + cov15[5][0] * f[0][2] + cov15[8][0] * f[0][3] + cov15[12][0]*f[0][4];
  w[3]=cov15[6][0] * f[0][0] + cov15[7][0] * f[0][1] + cov15[8][0] * f[0][2] + cov15[9][0] * f[0][3] + cov15[13][0]*f[0][4];
  w[4]=cov15[10][0]* f[0][0] + cov15[11][0]* f[0][1] + cov15[12][0]* f[0][2] + cov15[13][0]* f[0][3] + cov15[14][0]*f[0][4];

  covOut[0]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];

  w[0]=cov15[0][0] * f[1][0] + cov15[1][0] * f[1][1] + cov15[3][0] * f[1][2] + cov15[6][0] * f[1][3] + cov15[10][0]*f[1][4];
  w[1]=cov15[1][0] * f[1][0] + cov15[2][0] * f[1][1] + cov15[4][0] * f[1][2] + cov15[7][0] * f[1][3] + cov15[11][0]*f[1][4];
  w[2]=cov15[3][0] * f[1][0] + cov15[4][0] * f[1][1] + cov15[5][0] * f[1][2] + cov15[8][0] * f[1][3] + cov15[12][0]*f[1][4];
  w[3]=cov15[6][0] * f[1][0] + cov15[7][0] * f[1][1] + cov15[8][0] * f[1][2] + cov15[9][0] * f[1][3] + cov15[13][0]*f[1][4];
  w[4]=cov15[10][0]* f[1][0] + cov15[11][0]* f[1][1] + cov15[12][0]* f[1][2] + cov15[13][0]* f[1][3] + cov15[14][0]*f[1][4];

  covOut[1]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
  covOut[2]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];

  w[0]=cov15[0][0] * f[2][0] + cov15[1][0] * f[2][1] + cov15[3][0] * f[2][2] + cov15[6][0] * f[2][3] + cov15[10][0]*f[2][4];
  w[1]=cov15[1][0] * f[2][0] + cov15[2][0] * f[2][1] + cov15[4][0] * f[2][2] + cov15[7][0] * f[2][3] + cov15[11][0]*f[2][4];
  w[2]=cov15[3][0] * f[2][0] + cov15[4][0] * f[2][1] + cov15[5][0] * f[2][2] + cov15[8][0] * f[2][3] + cov15[12][0]*f[2][4];
  w[3]=cov15[6][0] * f[2][0] + cov15[7][0] * f[2][1] + cov15[8][0] * f[2][2] + cov15[9][0] * f[2][3] + cov15[13][0]*f[2][4];
  w[4]=cov15[10][0]* f[2][0] + cov15[11][0]* f[2][1] + cov15[12][0]* f[2][2] + cov15[13][0]* f[2][3] + cov15[14][0]*f[2][4];

  covOut[3]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
  covOut[4]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];
  covOut[5]=w[0]*f[2][0] + w[1]*f[2][1] + w[2]*f[2][2] + w[3]*f[2][3] + w[4]*f[2][4];

  w[0]=cov15[0][0] * f[3][0] + cov15[1][0] * f[3][1] + cov15[3][0] * f[3][2] + cov15[6][0] * f[3][3] + cov15[10][0]*f[3][4];
  w[1]=cov15[1][0] * f[3][0] + cov15[2][0] * f[3][1] + cov15[4][0] * f[3][2] + cov15[7][0] * f[3][3] + cov15[11][0]*f[3][4];
  w[2]=cov15[3][0] * f[3][0] + cov15[4][0] * f[3][1] + cov15[5][0] * f[3][2] + cov15[8][0] * f[3][3] + cov15[12][0]*f[3][4];
  w[3]=cov15[6][0] * f[3][0] + cov15[7][0] * f[3][1] + cov15[8][0] * f[3][2] + cov15[9][0] * f[3][3] + cov15[13][0]*f[3][4];
  w[4]=cov15[10][0]* f[3][0] + cov15[11][0]* f[3][1] + cov15[12][0]* f[3][2] + cov15[13][0]* f[3][3] + cov15[14][0]*f[3][4];

  covOut[6]=w[0]*f[0][0] + w[1]*f[0][1] + w[2]*f[0][2] + w[3]*f[0][3] + w[4]*f[0][4];
  covOut[7]=w[0]*f[1][0] + w[1]*f[1][1] + w[2]*f[1][2] + w[3]*f[1][3] + w[4]*f[1][4];
  covOut[8]=w[0]*f[2][0] + w[1]*f[2][1] + w[2]*f[2][2] + w[3]*f[2][3] + w[4]*f[2][4];
  covOut[9]=w[0]*f[3][0] + w[1]*f[3][1] + w[2]*f[3][2] + w[3]*f[3][3] + w[4]*f[3][4];

    // as dPinv/dpar = 0 if par != Pinv
  covOut[10]=cov15[10][0]*f[0][0] + cov15[11][0]*f[0][1] + cov15[12][0]*f[0][2] + cov15[13][0]*f[0][3] + cov15[14][0]*f[0][4];
  covOut[11]=cov15[10][0]*f[1][0] + cov15[11][0]*f[1][1] + cov15[12][0]*f[1][2] + cov15[13][0]*f[1][3] + cov15[14][0]*f[1][4];
  covOut[12]=cov15[10][0]*f[2][0] + cov15[11][0]*f[2][1] + cov15[12][0]*f[2][2] + cov15[13][0]*f[2][3] + cov15[14][0]*f[2][4];
  covOut[13]=cov15[10][0]*f[3][0] + cov15[11][0]*f[3][1] + cov15[12][0]*f[3][2] + cov15[13][0]*f[3][3] + cov15[14][0]*f[3][4];

  covOut[14]=cov15[14][0];


  if(covOut[0][0] < 0 || covOut[2][0] < 0 ||
     covOut[5][0] < 0 || covOut[9][0] < 0 ||
     covOut[14][0]< 0){
    GFException exc("RKTrackRepXY::Extrap() ==> output cov. matrix is wrong",__LINE__,__FILE__);
    std::vector< TMatrixT<double> > matrices;
    matrices.push_back(stateIn);
    matrices.push_back(cov15);
    matrices.push_back(stateOut);
    matrices.push_back(covOut);
    exc.setMatrices("stateIn cov15 stateOut cov15Out",matrices);
    exc.setFatal();
    throw exc;
  }

  return s;
}

ClassImp(RKTrackRepXY)
