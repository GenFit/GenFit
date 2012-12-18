#include "GFAbsRecoHit.h"
#include "SlTrackRep.h"
#include <cmath>

SlTrackRep::SlTrackRep()
  :GFAbsTrackRep(4),_backw(0)
{

}

SlTrackRep::SlTrackRep(const TMatrixT<double>& _state, 
                       const TMatrixT<double>& sigma)
  :GFAbsTrackRep(4),_backw(0)
{
  fState=_state;
  for(int i=0;i<sigma.GetNrows();++i){
    fCov[i][i]=sigma[i][0];
  }
  setReferencePlane(GFDetPlane(TVector3(0,0,0),TVector3(1,0,0),TVector3(0,1,0)));
}

SlTrackRep::SlTrackRep(const TVector3& pos, const TVector3& dir) :GFAbsTrackRep(4),_backw(0)
{
  GFDetPlane d(pos,dir);
  setReferencePlane(d);
  fState[0][0]=0.;
  fState[1][0]=0.;
  fState[2][0]=0.;
  fState[3][0]=0.;
  fCov[0][0]=1.;
  fCov[1][1]=1.;
  fCov[2][2]=1.;
  fCov[3][3]=1.;
}

SlTrackRep::SlTrackRep(const GFDetPlane& dp,
		       const TMatrixT<double>& _state, 
                       const TMatrixT<double>& sigma)
  :GFAbsTrackRep(4),_backw(0)
{
  fState=_state;
  for(int i=0;i<sigma.GetNrows();++i){
    fCov[i][i]=sigma[i][0];
  }
  setReferencePlane(dp);
}

SlTrackRep::SlTrackRep(const TMatrixT<double>& _state, 
                       const TMatrixT<double>& sigma,
		       const double z)
  :GFAbsTrackRep(4),_backw(0)
{

  fState=_state;
  for(int i=0;i<sigma.GetNrows();++i){
    fCov[i][i]=sigma[i][0];
  }
  setReferencePlane(GFDetPlane(TVector3(0,0,z),TVector3(1,0,0),TVector3(0,1,0)));
}


SlTrackRep::~SlTrackRep()
{
  
}
double SlTrackRep::extrapolate(const GFDetPlane& pl,
                               TMatrixT<double>& statePred){
  statePred.ResizeTo(fDimension,1);
  TVector3 o=pl.getO();
  TVector3 u=pl.getU();
  TVector3 v=pl.getV();
  TVector3 w=u.Cross(v);
  TVector3 ofrom=fRefPlane.getO();
  TVector3 ufrom=fRefPlane.getU();
  TVector3 vfrom=fRefPlane.getV();
  TVector3 wfrom=ufrom.Cross(vfrom);
  TVector3 p1=ofrom+fState[0][0]*ufrom+fState[1][0]*vfrom;

  TVector3 dir=(fState[2][0] * ufrom + fState[3][0] * vfrom + wfrom);//in global coordinates
  
  if(w*dir==0){
    std::cerr<<"track paralell to detector plane"<<std::endl
	     <<"extrapolation impossible, aborting"<<std::endl;
    throw;
  }
  
  double t= ((( w * o) - (w * p1)) / 
	     (w * dir));
  
  double dist=t*dir.Mag();
  
  TVector3 p2=p1+t*dir;
  
  double state0=(p2-o)*u;
  double state1=(p2-o)*v;
  double state2=(dir*u)/(dir*w);
  double state3=(dir*v)/(dir*w);

  statePred[0]=state0;
  statePred[1]=state1;
  statePred[2]=state2;
  statePred[3]=state3;
  return dist;
}
double SlTrackRep::extrapolate(const GFDetPlane& pl,
                               TMatrixT<double>& statePred,
                               TMatrixT<double>& covPred)
{
  //std::cout<<std::endl<<std::endl<<"extrapolate to plane with state and cov"<<std::endl<<std::endl;
  statePred.ResizeTo(fDimension,1);
  covPred.ResizeTo(fDimension,fDimension);
  TVector3 o=pl.getO();
  TVector3 u=pl.getU();
  TVector3 v=pl.getV();
  TVector3 w=u.Cross(v);
  TVector3 ofrom=fRefPlane.getO();
  TVector3 ufrom=fRefPlane.getU();
  TVector3 vfrom=fRefPlane.getV();
  TVector3 wfrom=ufrom.Cross(vfrom);
  /*
  std::cout<<"ufrom ";
  ufrom.Print();
  std::cout<<"vfrom ";
  vfrom.Print();
  std::cout<<"wfrom ";
  wfrom.Print();
  */
  TVector3 p1=ofrom+fState[0][0]*ufrom+fState[1][0]*vfrom;
  /*
  std::cout<<"p1 ";
  p1.Print();
  */
  TVector3 dir=(fState[2][0] * ufrom + fState[3][0] * vfrom + wfrom);//in global coordinates
  /*
    std::cout<<"dir ";
    dir.Print();
  */
  if(w*dir==0){
    std::cerr<<"track paralell to detector plane"<<std::endl
	     <<"extrapolation impossible, aborting"<<std::endl;
    throw;
  }

  double t= ((( w * o) - (w * p1)) / 
	     (w * dir));
  
  double dist=t*dir.Mag();
  
  TVector3 p2=p1+t*dir;
  /*
    std::cout<<"xtrapolate point p2 ";
    p2.Print();
  */
  double state0=(p2-o)*u;
  double state1=(p2-o)*v;
  double state2=(dir*u)/(dir*w);
  double state3=(dir*v)/(dir*w);
  /*
  std::cout<<"start state:"<<std::endl;
  std::cout<<fState[0][0]<<std::endl
	   <<fState[1][0]<<std::endl
	   <<fState[2][0]<<std::endl
	   <<fState[3][0]<<std::endl;
  std::cout<<"end state:"<<std::endl;
  std::cout<<fState0<<std::endl
	   <<fState1<<std::endl
	   <<fState2<<std::endl
	   <<fState3<<std::endl<<std::endl;
  */
  TMatrixT<double> jacobian(4,4);
  /*
  cg3 = .(ufrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ufrom) * c * ufrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ufrom) * d * vfrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ufrom) *wfrom, u);

  cg4 = .(vfrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, vfrom) * c * ufrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, vfrom) * d * vfrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, vfrom) *wfrom, u);

  cg5 = .(- pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * c * ufrom * .(w, ufrom) 
	  + 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, o) * ufrom 
	  - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * d * vfrom * .(w, ufrom) 
	  - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) *wfrom * .(w, ufrom) 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * c * ufrom * .(w, ufrom) 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ofrom + a * ufrom + b * vfrom) * ufrom 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * d * vfrom * .(w, ufrom) 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) *wfrom * .(w, ufrom), u);

  cg6 = .(-pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * c * ufrom * .(w, vfrom) 
	  - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * d * vfrom * .(w, vfrom) 
	  + 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, o) * vfrom 
	  - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) *wfrom * .(w, vfrom) 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * c * ufrom * .(w, vfrom) 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * d * vfrom * .(w, vfrom) 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ofrom + a * ufrom + b * vfrom) * vfrom 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) *wfrom * .(w, vfrom), u);

  cg7 = .(ufrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ufrom) * c * ufrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ufrom) * d * vfrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ufrom) *wfrom, v);

  cg8 = .(vfrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, vfrom) * c * ufrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, vfrom) * d * vfrom 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, vfrom) *wfrom, v);

  cg9 = .(-pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * c * ufrom * .(w, ufrom) 
	  + 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, o) * ufrom 
	  - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * d * vfrom * .(w, ufrom) 
	  - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) *wfrom * .(w, ufrom) 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * c * ufrom * .(w, ufrom) 
	  - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ofrom + a * ufrom + b * vfrom) * ufrom 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * d * vfrom * .(w, ufrom) 
	  + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) *wfrom * .(w, ufrom), v);

  cg10 = .(-pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * c * ufrom * .(w, vfrom) 
	   - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) * d * vfrom * .(w, vfrom) 
	   + 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, o) * vfrom 
	   - pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, o) *wfrom * .(w, vfrom)
	   + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * c * ufrom * .(w, vfrom) 
	   + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) * d * vfrom * .(w, vfrom) 
	   - 1 / .(w, c * ufrom + d * vfrom +wfrom) * .(w, ofrom + a * ufrom + b * vfrom) * vfrom 
	   + pow(.(w, c * ufrom + d * vfrom +wfrom), -2) * .(w, ofrom + a * ufrom + b * vfrom) *wfrom * .(w, vfrom), v);
 
  cg11 = 0;
  cg12 = 0;
 
  cg13 = .(ufrom, u) / .(c * ufrom + d * vfrom +wfrom, w) 
    - .(c * ufrom + d * vfrom +wfrom, u) * pow(.(c * ufrom + d * vfrom +wfrom, w), -2) * .(ufrom, w);

  cg14 = .(vfrom, u) / .(c * ufrom + d * vfrom +wfrom, w) 
    - .(c * ufrom + d * vfrom +wfrom, u) * pow(.(c * ufrom + d * vfrom +wfrom, w), -2) * .(vfrom, w);
  cg15 = 0;
  cg16 = 0;

  cg17 = .(ufrom, v) / .(c * ufrom + d * vfrom +wfrom, w) 
    - .(c * ufrom + d * vfrom +wfrom, v) * pow(.(c * ufrom + d * vfrom +wfrom, w), -2) * .(ufrom, w);
    
    cg18 = .(vfrom, v) / .(c * ufrom + d * vfrom +wfrom, w) 
    - .(c * ufrom + d * vfrom +wfrom, v) * pow(.(c * ufrom + d * vfrom +wfrom, w), -2) * .(vfrom, w);

  */


  double jacobian0 = (ufrom 
		      - 1 / (w * dir) * (w * ufrom) * fState[2][0] * ufrom 
		      - 1 / (w * dir) * (w * ufrom) * fState[3][0] * vfrom 
		      - 1 / (w * dir) * (w * ufrom) * wfrom) * u;

  double jacobian1 = (vfrom 
		      - 1 / (w * dir) * (w * vfrom) * fState[2][0] * ufrom 
		      - 1 / (w * dir) * (w * vfrom) * fState[3][0] * vfrom 
		      - 1 / (w * dir) * (w * vfrom) * wfrom)* u;

  double jacobian2 = (- pow(  (w * dir), -2) * (w * o) * fState[2][0] * ufrom * (w * ufrom) 
		      + 1 / (w * dir) * (w * o) * ufrom 
		      - pow((w * dir), -2) * (w * o) * fState[3][0] * vfrom * (w * ufrom) 
		      - pow((w * dir), -2) * (w * o) * wfrom * (w * ufrom) 
		      + pow((w * dir), -2) * (w * p1) * fState[2][0] * ufrom * (w * ufrom) 
		      - 1 / (w * dir) * (w * p1) * ufrom 
		      + pow((w * dir), -2) * (w * p1) * fState[3][0] * vfrom * (w * ufrom) 
		      + pow((w * dir), -2) * (w * p1) * wfrom * (w * ufrom))* u;

  double jacobian3 = (-pow( ( w * dir), -2) * (w * o) * fState[2][0] * ufrom * (w * vfrom) 
		      - pow((w * dir), -2) * (w * o) * fState[3][0] * vfrom * (w * vfrom) 
		      + 1 / (w * dir) * (w * o) * vfrom 
		      - pow((w * dir), -2) * (w * o) * wfrom * (w * vfrom) 
		      + pow((w * dir), -2) * (w * p1) * fState[2][0] * ufrom * (w * vfrom) 
		      + pow((w * dir), -2) * (w * p1) * fState[3][0] * vfrom * (w * vfrom) 
		      - 1 / (w * dir) * (w * p1) * vfrom 
		      + pow((w * dir), -2) * (w * p1) * wfrom * (w * vfrom))* u;

  double jacobian4 = (ufrom 
		      - 1 / (w * dir) * (w * ufrom) * fState[2][0] * ufrom 
		      - 1 / (w * dir) * (w * ufrom) * fState[3][0] * vfrom 
		      - 1 / (w * dir) * (w * ufrom) * wfrom)* v;

  double jacobian5 = (vfrom 
		      - 1 / (w * dir) * (w * vfrom) * fState[2][0] * ufrom 
		      - 1 / (w * dir) * (w * vfrom) * fState[3][0] * vfrom 
		      - 1 / (w * dir) * (w * vfrom) * wfrom)* v;

  double jacobian6 = (- pow( ( w * dir), -2) * (w * o) * fState[2][0] * ufrom * (w * ufrom) 
		      + 1 / (w * dir) * (w * o) * ufrom 
		      - pow((w * dir), -2) * (w * o) * fState[3][0] * vfrom * (w * ufrom) 
		      - pow((w * dir), -2) * (w * o) * wfrom * (w * ufrom) 
		      + pow((w * dir), -2) * (w * p1) * fState[2][0] * ufrom * (w * ufrom) 
		      - 1 / (w * dir) * (w * p1) * ufrom 
		      + pow((w * dir), -2) * (w * p1) * fState[3][0] * vfrom * (w * ufrom) 
		      + pow((w * dir), -2) * (w * p1) * wfrom * (w * ufrom))* v;

  double jacobian7 = (- pow( ( w * dir), -2) * (w * o) * fState[2][0] * ufrom * (w * vfrom) 
		      - pow((w * dir), -2) * (w * o) * fState[3][0] * vfrom * (w * vfrom) 
		      + 1 / (w * dir) * (w * o) * vfrom 
		      - pow((w * dir), -2) * (w * o) * wfrom * (w * vfrom) 
		      + pow((w * dir), -2) * (w * p1) * fState[2][0] * ufrom * (w * vfrom) 
		      + pow((w * dir), -2) * (w * p1) * fState[3][0] * vfrom * (w * vfrom) 
		      - 1 / (w * dir) * (w * p1) * vfrom 
		      + pow((w * dir), -2) * (w * p1) * wfrom * (w * vfrom))* v;
  double jacobian8 = 0;
  double jacobian9 = 0;
  
  double jacobian10 = ((ufrom * u) / (w * dir) 
		       - (dir * u) * pow((w * dir), -2) * (w * ufrom));
  
  double jacobian11 = ((vfrom * u) / (w * dir) 
		       - (dir * u) * pow((w * dir), -2) * (w * vfrom));
  double jacobian12 = 0;
  double jacobian13 = 0;
  
  double jacobian14 = ((ufrom * v) / (w * dir) 
		       - (dir * v) * pow((w * dir), -2) * (w * ufrom));
  
  double jacobian15 = ((vfrom * v) / (w * dir) 
		       - (dir * v) * pow((w * dir), -2) * (w * vfrom));
 
 jacobian[0][0]=jacobian0; jacobian[0][1]=jacobian1; jacobian[0][2]=jacobian2; jacobian[0][3]=jacobian3;
 jacobian[1][0]=jacobian4; jacobian[1][1]=jacobian5; jacobian[1][2]=jacobian6; jacobian[1][3]=jacobian7;
 jacobian[2][0]=jacobian8; jacobian[2][1]=jacobian9; jacobian[2][2]=jacobian10; jacobian[2][3]=jacobian11;
 jacobian[3][0]=jacobian12; jacobian[3][1]=jacobian13; jacobian[3][2]=jacobian14; jacobian[3][3]=jacobian15;
 TMatrixT<double> jacobianT(jacobian);
 jacobianT.T();
 covPred=jacobian*fCov*(jacobianT);
 // std::cout<<"covariance[0][0]"<<covPred[0][0]<<std::endl;
 statePred[0]=state0;
 statePred[1]=state1;
 statePred[2]=state2;
 statePred[3]=state3;
 return dist;
}
double SlTrackRep::extrapolateToPoint(const TVector3& pos,
				     TVector3& poca,
				     TVector3& dirInPoca){

  TVector3 ofrom=fRefPlane.getO();
  TVector3 ufrom=fRefPlane.getU();
  TVector3 vfrom=fRefPlane.getV();
  TVector3 wfrom=ufrom.Cross(vfrom);
  TVector3 pfrom=ofrom + fState[0][0]*ufrom + fState[1][0]*vfrom;
  TVector3 dir=(fState[2][0] * ufrom + fState[3][0] * vfrom + wfrom);//in global coordinates
  //  pfrom.Print();
  dir=dir.Unit();
  //dir.Print();
  double t = (dir * pos - dir * pfrom) / (dir * dir);
  
  poca=pfrom + t * dir;
  dirInPoca=dir.Unit();
  return t;
}
double SlTrackRep::extrapolateToLine(const TVector3& point1,
				   const TVector3& point2,
				   TVector3& poca,
				   TVector3& dirInPoca,
				   TVector3& poca_onwire){
  
  TVector3 ofrom=fRefPlane.getO();
  TVector3 ufrom=fRefPlane.getU();
  TVector3 vfrom=fRefPlane.getV();
  TVector3 wfrom=ufrom.Cross(vfrom);

  TVector3 pfrom=ofrom + fState[0][0]*ufrom + fState[1][0]*vfrom;
  TVector3 dir=(fState[2][0] * ufrom + fState[3][0] * vfrom + wfrom);//in global coordinates
  TVector3 lineDir=point2-point1;
  //normal normal to both lineDir and Dir
  TVector3 normal(dir.y()*lineDir.z()-dir.z()*lineDir.y(),
		  dir.x()*lineDir.z()-dir.z()*lineDir.x(),
		  dir.x()*lineDir.y()-dir.y()*lineDir.x());
  
  normal=normal.Unit();
  TVector3 planeNorm=dir.Cross(normal);
  double t=(planeNorm*point2-planeNorm*pfrom)/(planeNorm*dir);
  poca=pfrom+t*dir;
  double t2=(lineDir*poca - lineDir*point1)/(lineDir*lineDir);
  poca_onwire=point1+lineDir*t2;
  dirInPoca=dir;
  return t;
}
TVector3 SlTrackRep::getPos(const GFDetPlane& pl){
  TMatrixT<double> statePred(fState);
  if(pl!=fRefPlane){
    extrapolate(pl,statePred);
  }
  return pl.getO()+(statePred[0][0]*pl.getU())+(statePred[1][0]*pl.getV());
}
TVector3 SlTrackRep::getMom(const GFDetPlane& pl){
  TMatrixT<double> statePred(fState);
  if(pl!=fRefPlane){
    extrapolate(pl,statePred);
  }
  TVector3 ret = pl.getNormal()+(statePred[2][0]*pl.getU())+(statePred[3][0]*pl.getV());
  ret.SetMag(1.);
  return ret;
}
void SlTrackRep::getPosMom(const GFDetPlane& pl,TVector3& pos,TVector3& mom){
  pos=getPos(pl);
  mom=getMom(pl);
}
