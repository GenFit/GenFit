//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class GeaneTrackRep
//      see GeaneTrackRep.hh for details
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "GeaneTrackRep2.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <cmath>
#include <cassert>

// Collaborating Class Headers --------
#include "GFAbsRecoHit.h"
#include "GFException.h"
#include "TGeant3.h"
#include "TDatabasePDG.h"
// Class Member definitions -----------

#define THETACUT 0.4

GeaneTrackRep2::GeaneTrackRep2()
  : GFAbsTrackRep(6), fPdg(0)
{

}

GeaneTrackRep2::GeaneTrackRep2(const GFDetPlane& plane,
			       const TVector3& mom,
			       const TVector3& poserr,
			       const TVector3& momerr,
			       int PDGCode) 
  : GFAbsTrackRep(6), fPdg(PDGCode)
{

  fG3ParticleID=TDatabasePDG::Instance()->ConvertPdgToGeant3(fPdg);
  //TDatabasePDG has charge in units of |e|/3
  fCharge=TDatabasePDG::Instance()->GetParticle(fPdg)->Charge()/3.;
  /*
  std::cout << "Initializing GeaneTrackRep2 for particle "
	    << TDatabasePDG::Instance()->GetParticle(fPdg)->GetName()
	    << std::endl;
  */
  fRefPlane=plane;
  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);

  fState[3][0] = 0.;
  fState[4][0] = 0.;
  fState[0][0] = fCharge/mom.Mag();
  double pu,pv,pw;
  pu = mom*u;
  pv = mom*v;
  pw = mom*w;

  if(fabs(pw)<1.e-10){
    GFException exc("fabs(pw<1.e-10)",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  fState[1][0] = pu/pw;
  fState[2][0] = pv/pw;
  if(pw>0.){
    //set spu
    fState[5][0] = 1.;
  }
  else{
    //set spu
    fState[5][0] = -1.;
  }
  fCov[3][3] = poserr.X()*poserr.X() * u.X()*u.X() +
    poserr.Y()*poserr.Y() * u.Y()*u.Y() +
    poserr.Z()*poserr.Z() * u.Z()*u.Z();
  fCov[4][4] = poserr.X()*poserr.X() * v.X()*v.X() +
    poserr.Y()*poserr.Y() * v.Y()*v.Y() +
    poserr.Z()*poserr.Z() * v.Z()*v.Z();
  fCov[0][0] = fCharge*fCharge/pow(mom.Mag(),6.) * 
    (mom.X()*mom.X() * momerr.X()*momerr.X()+
     mom.Y()*mom.Y() * momerr.Y()*momerr.Y()+
     mom.Z()*mom.Z() * momerr.Z()*momerr.Z());
  fCov[1][1] = pow((u.X()/pw - w.X()*pu/(pw*pw)),2.)*momerr.X()*momerr.X() +
    pow((u.Y()/pw - w.Y()*pu/(pw*pw)),2.)*momerr.Y()*momerr.Y() +
    pow((u.Z()/pw - w.Z()*pu/(pw*pw)),2.)*momerr.Z()*momerr.Z();
  fCov[2][2] = pow((v.X()/pw - w.X()*pv/(pw*pw)),2.)*momerr.X()*momerr.X() +
    pow((v.Y()/pw - w.Y()*pv/(pw*pw)),2.)*momerr.Y()*momerr.Y() +
    pow((v.Z()/pw - w.Z()*pv/(pw*pw)),2.)*momerr.Z()*momerr.Z();
  //fCov[5][5] is zero - spu has no error

}


GeaneTrackRep2::~GeaneTrackRep2()
{
  
}




double
GeaneTrackRep2::extrapolate(const GFDetPlane& pl, 
			   TMatrixT<double>& statePred)
{
  TMatrixT<double> covPred(6,6);
  return  extrapolate(pl,statePred,covPred);
  //! TODO: make this faster by neglecting covariances ?
}


double
GeaneTrackRep2::extrapolate(const GFDetPlane& pl, 
			   TMatrixT<double>& statePred,
			   TMatrixT<double>& covPred)
{
  if(fabs(getMom(fRefPlane).Theta()/TMath::Pi()*180.) < THETACUT){
    GFException exc("GEANE propagation not possible for p.theta<THETACUT",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  TGeant3 *gMC3 = (TGeant3*) gMC;        
  if(gMC3==NULL){
    std::cerr << "GeaneTrackRep2: TGeant3 has not been initialized. -> abort"
	      << std::endl;
    throw;
  }

  statePred.ResizeTo(fDimension,1);covPred.ResizeTo(fDimension,fDimension);
  TVector3 oto=pl.getO();TVector3 uto=pl.getU();TVector3 vto=pl.getV();

  TVector3 ofrom=fRefPlane.getO();TVector3 ufrom=fRefPlane.getU();TVector3 vfrom=fRefPlane.getV();
  TVector3 wfrom=ufrom.Cross(vfrom);;

  Float_t pli[6];
  pli[0]=ufrom.X();pli[1]=ufrom.Y();pli[2]=ufrom.Z();pli[3]=vfrom.X();pli[4]=vfrom.Y();pli[5]=vfrom.Z();

  Float_t plf[12];
  plf[0]=uto.X();plf[1]=uto.Y();plf[2]=uto.Z();plf[3]=vto.X();plf[4]=vto.Y();plf[5]=vto.Z();plf[6]=oto.X();plf[7]=oto.Y();plf[8]=oto.Z();

  TString geaneOption;
  geaneOption="PE";

  Float_t ein[15];
  int count=0;;
  for(int i=0; i<5;++i){
    for(int j=i;j<5;++j){
      ein[count++]=fCov[i][j];
    }
  }
  //convert covariance from q/p to 1/p
  ein[0] /= fCharge*fCharge;ein[1] *= fCharge;ein[2] *= fCharge;ein[3] *= fCharge;ein[4] *= fCharge;

  if(fState[3][0]==0)fState[3][0]=1E-4;
  if(fState[4][0]==0)fState[4][0]=1E-4;
  gMC3->Eufilp(1, ein, pli, plf);

  TVector3 x1vect;
  Float_t x1[3];Float_t x2[3];Float_t p1[3];Float_t p2[3];

  x1vect = ofrom + fState[3][0]*ufrom + fState[4][0]*vfrom;
  x1vect.GetXYZ(x1);
  double mom = 1.e30;
  //double pOVERpw = sqrt(1.+fState[1][0]*fState[1][0]+fState[2][0]*fState[2][0]);
  if(fabs(fState[0][0])>1.e-30) mom = 1./fabs(fState[0][0]);

  //momentum in 'local mars'
  TVector3 p1vect = fState[5][0]*(wfrom+fState[1][0]*ufrom+fState[2][0]*vfrom);
  p1vect.SetMag(mom);

  p1vect.GetXYZ(p1);

  //Determine whether we are doing a forward or backward step:
  TVector3 dir=pl.dist(x1vect); // direction from pos to plane;
  if((dir*p1vect)<0){
    geaneOption="BPE";
  }

  for(unsigned int i=0;i<3;++i) {
    x2[i]=0.;
    p2[i]=0.;
  }

  gMC3->Ertrak(x1,p1,x2,p2,fG3ParticleID, geaneOption.Data());

  if(x2[0]<-1.E29) {
    GFException e("obsolete exception?: x2[0]<-1.E29",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  if(fabs(p2[0])==0. && fabs(p2[1])==0. && fabs(p2[2])==0) {
    GFException e("fabs(p2[0])==0. && fabs(p2[1])==0. && fabs(p2[2])==0",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  double trklength = gMC3->TrackLength();
  Ertrio_t *ertrio = gMC3->fErtrio;

#ifndef OLDVMC
  Ertrio1_t *ertrio1 = gMC3->fErtrio1;
  if(ertrio1->iertr != 0){
    GFException e("GEANE error flag for numerical divergence detected",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
#endif

  Double_t cov15[15];

  for(Int_t i=0;i<15;i++) {
    cov15[i]=ertrio->errout[i]; 
    //go from 1/p to q/p
    if(i == 0) cov15[i] = cov15[i] * fCharge * fCharge;
    if(i > 0 && i < 5) cov15[i] = cov15[i] * fCharge;
  }

  count=0;
  for(int i=0;i<5;++i){
    for(int j=i;j<5;++j){
      covPred[i][j]=cov15[count];
      if(i!=j)covPred[j][i]=cov15[count];
      ++count;
    }
  }
  TVector3 d(x2[0]-oto.X(),x2[1]-oto.Y(),x2[2]-oto.Z());
  statePred[3][0] = d*uto;
  statePred[4][0] = d*vto;
  TVector3 p2vect(p2[0],p2[1],p2[2]);
  statePred[0][0] = fCharge/p2vect.Mag();
  double pu,pv,pw;
  pu = p2vect*uto;
  pv = p2vect*vto;
  pw = p2vect*(uto.Cross(vto));

  if(fabs(pw)<1.e-10){
    GFException exc("fabs(pw)<1.e-10",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  if(pw>0.){
    statePred[5][0]=1.;
    statePred[1][0] = pu/fabs(pw);
    statePred[2][0] = pv/fabs(pw);
  }
  else{
    statePred[5][0]=-1.;
    statePred[1][0] = -1.*pu/fabs(pw);
    statePred[2][0] = -1.*pv/fabs(pw);
  }
  return trklength;
}

void
GeaneTrackRep2::extrapolateToPoint(const TVector3& pos,
				 TVector3& poca,
				 TVector3& normVec){
  
  if(fabs(getMom(fRefPlane).Theta()/TMath::Pi()*180.) < THETACUT){
    GFException exc("GEANE propagation not possible for p.theta<THETACUT",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  TGeant3 *gMC3 = (TGeant3*) gMC;        
  if(gMC3==NULL){
    std::cerr << "GeaneTrackRep2: TGeant3 has not been initialized. -> abort"
	      << std::endl;
    throw;
  }


  Float_t nullVec[3];
  nullVec[0]=0.;  nullVec[1]=0.;  nullVec[2]=0.;
  Float_t posArr[3];
  pos.GetXYZ(posArr);
  //  gMC3->SetClose(pca,pf,dst,w1,w2,po1,po2,po3,cl);
  gMC3->SetClose(1,posArr,999,nullVec,nullVec,nullVec,nullVec,nullVec,nullVec);

  
  Float_t ein[15];
  int count=0;;
  for(int i=0; i<5;++i){
    for(int j=i;j<5;++j){
      ein[count++]=fCov[i][j];
    }
  }
  //convert covariance from q/p to 1/p
  ein[0] /= fCharge*fCharge;ein[1] *= fCharge;ein[2] *= fCharge;ein[3] *= fCharge;ein[4] *= fCharge;
  
  TVector3 mypos = getPos(fRefPlane);
  TVector3 mymom = getMom(fRefPlane);
  Float_t x1[3];
  Float_t p1[3];
  Float_t x2[3];
  Float_t p2[3];
  Float_t po1[3];
  Float_t po2[3];
  Float_t po3[3];
  Float_t clen[3];
  for(unsigned int i=0;i<3;++i){
    x2[i]=0.;p2[i]=0.;po1[i]=0.;po2[i]=0.;po3[i]=0.;clen[i]=0.;
  }
  mypos.GetXYZ(x1);
  mymom.GetXYZ(p1);
  Float_t maxlen[1] = {2.*((mypos-pos).Mag())};
  gMC3->Eufill(1, ein, maxlen);


  TString geaneOption;
  geaneOption="LO";
  //Determine whether we are doing a forward or backward step:
  TVector3 dir=fRefPlane.dist(pos); // direction from pos to plane;
  if((dir*mymom)>0){
    geaneOption="BLO";
  }

  gMC3->Ertrak(x1,p1,x2,p2,fG3ParticleID, geaneOption.Data());

  gMC3->GetClose(po1,po2,po3,clen);

  TVector3 poV1(po1);
  TVector3 poV2(po2);
  TVector3 poV3(po3);
  if((fabs(clen[0])<1.E-8 && fabs(clen[1])<1.E-8) ||
     fabs(clen[0]-clen[1])<1.E-8){
    poca2Line(poV2,poV3,pos,poca);
    normVec = poV3-poV2;
  }

  else if(fabs(clen[1]-clen[2])<1.E-8){
    poca2Line(poV1,poV2,pos,poca);
    normVec = poV2-poV1;
  }
  else{
    TVector3 result12,result23;
    poca2Line(poV1,poV2,pos,result12);
    poca2Line(poV2,poV3,pos,result23);
    if((result12-pos).Mag()<(result23-pos).Mag()){
      poca = result12;
      if( (poV2-poca).Mag() > 0.25*(poV1-poV2).Mag() ){
	normVec = poV2-poV1;
      }
      else{//poca is close to poV2, so take poV3-poV1 as direction vector
	normVec = poV3-poV1;
      }
    }
    else{
      poca = result23;
      if( (poV2-poca).Mag() > 0.25*(poV2-poV3).Mag() ){
	normVec = poV3-poV2;
      }
      else{//poca is close to poV2, so take poV3-poV1 as direction vector
	normVec = poV3-poV1;
      }
    }
  }
  normVec.SetMag(1.);
  /*
  int dim = getDim();
  TMatrixT<double> statePred(dim,1);
  TMatrixT<double> covPred(dim,dim);
  //std::cout<<"GeaneTrackRep2::extrapolateToPoint"<<std::endl;
  //fRefPlane.Print();

  TVector3 ofrom=fRefPlane.getO();
  TVector3 ufrom=fRefPlane.getU();
  TVector3 vfrom=fRefPlane.getV();

  _geane->SetPoint(pos);
  _geane->PropagateFromPlane(ufrom,vfrom);

  double cova[15];
  int count=0;;
  for(int i=0; i<5;++i){
    for(int j=i;j<5;++j){
      cova[count++]=cov[i][j];
    }
  }
  // protect against low momentum:
  if(fabs(state[0][0])>10){
    GFException exc("GeaneTrackRep2: PROTECT AGAINST LOW MOMENTA",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // protect against (x,y)=(0,0)
  if(state[3][0]==0)state[3][0]=1E-4;
  if(state[4][0]==0)state[4][0]=1E-4;
  
  FairTrackParP par(state[3][0],state[4][0],state[1][0],state[2][0],state[0][0],cova,ofrom,ufrom,vfrom,_spu);
  //par.Print();
  bool backprop=_backw<0;
  if(_backw==0){
    // check if new point is after or before my position
    TVector3 dir=fRefPlane.dist(pos); // direction from pos to plane;
	backprop= (dir*getMom(fRefPlane))>0;
  }
  if(!backprop){ // point lies in same direction of flight as momentum
    //std::cout<<" Propagate in flight direction"<<std::endl;
    _geane->PropagateToVirtualPlaneAtPCA(1);
  }
  else{
    //std::cout<<" backPropagate"<<std::endl;
    _geane->BackTrackToVirtualPlaneAtPCA(1);
  }

  FairTrackParP result;
  Bool_t prop = kTRUE;
  prop = _geane->Propagate(&par,&result,fPdg);   //211
  if (prop==kFALSE) {
    GFException exc("GEANE propagation failed",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
    //pl=fRefPlane;
    //return pos;
  }

  statePred[0][0]=result.GetQp();
  statePred[1][0]=result.GetTV();
  statePred[2][0]=result.GetTW();
  statePred[3][0]=result.GetV();
  statePred[4][0]=result.GetW();

  double* rescov=result.GetCov();
  count=0;
  for(int i=0;i<5;++i){
    for(int j=i;j<5;++j){
      covPred[i][j]=rescov[count];
      if(i!=j)covPred[j][i]=rescov[count];
      ++count;
    }
  }

  poca.SetXYZ(result.GetX(),result.GetY(),result.GetZ());
  normVec = result.GetJVer().Cross( result.GetKVer() );
  */
}


void 
GeaneTrackRep2::extrapolateToLine(const TVector3& point1,
				 const TVector3& point2,
				 TVector3& poca,
				 TVector3& normVec,
				 TVector3& poca_onWire)
{
  if(fabs(getMom(fRefPlane).Theta()/TMath::Pi()*180.) < THETACUT){
    GFException exc("GEANE propagation not possible for p.theta<THETACUT",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  TGeant3 *gMC3 = (TGeant3*) gMC;        
  if(gMC3==NULL){
    std::cerr << "GeaneTrackRep2: TGeant3 has not been initialized. -> abort"
	      << std::endl;
    throw;
  }


  Float_t nullVec[3];
  nullVec[0]=0.;  nullVec[1]=0.;  nullVec[2]=0.;
  Float_t w1Arr[3];
  Float_t w2Arr[3];
  point1.GetXYZ(w1Arr);
  point2.GetXYZ(w2Arr);
  //  gMC3->SetClose(pca,pf,dst,w1,w2,po1,po2,po3,cl);
  gMC3->SetClose(2,nullVec,999,w1Arr,w2Arr,nullVec,nullVec,nullVec,nullVec);

  
  Float_t ein[15];
  int count=0;;
  for(int i=0; i<5;++i){
    for(int j=i;j<5;++j){
      ein[count++]=fCov[i][j];
    }
  }
  //convert covariance from q/p to 1/p
  ein[0] /= fCharge*fCharge;ein[1] *= fCharge;ein[2] *= fCharge;ein[3] *= fCharge;ein[4] *= fCharge;
  
  TVector3 mypos = getPos(fRefPlane);
  TVector3 mymom = getMom(fRefPlane);
  Float_t x1[3];
  Float_t p1[3];
  Float_t x2[3];
  Float_t p2[3];
  Float_t po1[3];
  Float_t po2[3];
  Float_t po3[3];
  Float_t clen[3];
  for(unsigned int i=0;i<3;++i){
    x2[i]=0.;p2[i]=0.;po1[i]=0.;po2[i]=0.;po3[i]=0.;clen[i]=0.;
  }
  mypos.GetXYZ(x1);
  mymom.GetXYZ(p1);

  TVector3 pointOnWireClosestToMyPos;
  poca2Line(point1,point2,mypos,pointOnWireClosestToMyPos);

  Float_t maxlen[1] = {2.*((mypos-pointOnWireClosestToMyPos).Mag())};
  gMC3->Eufill(1, ein, maxlen);

  //point1 muss pocaOnWire sein
  TString geaneOption;
  geaneOption="LO";
  //Determine whether we are doing a forward or backward step:
  TVector3 dir=fRefPlane.dist(pointOnWireClosestToMyPos); // direction from pos to plane;
  if((dir*mymom)>0){
    geaneOption="BLO";
  }

  gMC3->Ertrak(x1,p1,x2,p2,fG3ParticleID, geaneOption.Data());

  gMC3->GetClose(po1,po2,po3,clen);

  TVector3 poV1(po1);
  TVector3 poV2(po2);
  TVector3 poV3(po3);
  if((fabs(clen[0])<1.E-8 && fabs(clen[1])<1.E-8) ||
     fabs(clen[0]-clen[1])<1.E-8){
    pocaLine2Line(poV2,poV3-poV2,point1,point2-point1,poca,poca_onWire);
    normVec = poV3-poV2;
  }
  else  if(fabs(clen[1]-clen[2])<1.E-8){
    pocaLine2Line(poV1,poV2-poV1,point1,point2-point1,poca,poca_onWire);
    normVec = poV2-poV1;
  }
  else{
    TVector3 result12_1,result12_2,result23_1,result23_2;
    pocaLine2Line(poV1,poV2-poV1,point1,point2-point1,result12_1,result12_2);
    pocaLine2Line(poV2,poV3-poV2,point1,point2-point1,result23_1,result23_2);
    if((result12_1-result12_2).Mag()<(result23_1-result23_2).Mag()){
      poca = result12_1;
      poca_onWire = result12_2;
      if( (poV2-poca).Mag() > 0.25*(poV1-poV2).Mag() ){
	normVec = poV2-poV1;
      }
      else{//poca is close to poV2, so take poV3-poV1 as direction vector
	normVec = poV3-poV1;
      }
    }
    else{
      poca = result23_1;
      poca_onWire = result23_2;
      if( (poV2-poca).Mag() > 0.25*(poV2-poV3).Mag() ){
	normVec = poV3-poV2;
      }
      else{//poca is close to poV2, so take poV3-poV1 as direction vector
	normVec = poV3-poV1;
      }
    }
  }
  normVec.SetMag(1.);


}

void GeaneTrackRep2::poca2Line(const TVector3& extr1,const TVector3& extr2,const TVector3& point,TVector3& result){
  
  TVector3 theWire = extr2-extr1;
  if(theWire.Mag()<1.E-8){
    GFException exc("GeaneTrackRep2::poca2Line(): try to find poca between line and point, but the line is really just a point",__LINE__,__FILE__);
    throw exc;
  }
  double t = 1./(theWire*theWire)*(point*theWire+extr1*extr1-extr1*extr2);
  result = extr1+t*theWire;
}

void GeaneTrackRep2::pocaLine2Line(const TVector3& point1,const TVector3& line1,const TVector3& point2, const TVector3& line2,TVector3& result1,TVector3& result2){
  TVector3 l(line1);
  l.SetMag(1.);
  TVector3 k(line2);
  k.SetMag(1.);
  double kl = k*l;

  if(fabs(kl)<1.E-8){
    GFException exc("GeaneTrackRep2::pocaLine2Line(): lines are parallel",__LINE__,__FILE__);
    throw exc;
  }
  double s = 1./(kl*kl-1)*(point2*k-point1*k-(point2*l-point1*l)*kl);
  double t = (point2*l+s*kl-point1*l);
  result1 = point1+t*l;
  result2 = point2+s*k;
}



TVector3 
GeaneTrackRep2::getPos(const GFDetPlane& pl)
{
  TMatrixT<double> statePred(fState);
  if(pl!=fRefPlane)extrapolate(pl,statePred);
  return pl.getO()+(statePred[3][0]*pl.getU())+(statePred[4][0]*pl.getV());
}
 
TVector3 
GeaneTrackRep2::getMom(const GFDetPlane& pl)
{
  TMatrixT<double> statePred(fState);
  if(pl!=fRefPlane)extrapolate(pl,statePred);
  //pl.Print();
  //statePred.Print();

  TVector3 mom = statePred[5][0]*(pl.getNormal()+statePred[1][0]*pl.getU()+statePred[2][0]*pl.getV());
  mom.SetMag(1./fabs(statePred[0][0]));
  return mom;
}
void
GeaneTrackRep2::getPosMom(const GFDetPlane& pl,TVector3& pos, TVector3& mom)
{
  TMatrixT<double> statePred(fState);
  if(pl!=fRefPlane)extrapolate(pl,statePred);
  mom = statePred[5][0]*(pl.getNormal()+statePred[1][0]*pl.getU()+statePred[2][0]*pl.getV());
  mom.SetMag(1./fabs(statePred[0][0]));
  pos = pl.getO()+(statePred[3][0]*pl.getU())+(statePred[4][0]*pl.getV());
}


void
GeaneTrackRep2::getPosMomCov(const GFDetPlane& pl,TVector3& pos,TVector3& mom,TMatrixT<double>& cov){
  cov.ResizeTo(6,6);
  std::cerr<<"insert brain here " << __FILE__ << " " << __LINE__
	   << " ->abort" <<std::endl;
  throw;
}
 

ClassImp(GeaneTrackRep2)
