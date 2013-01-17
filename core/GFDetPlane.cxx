/* Copyright 2008-2010, Technische Universitaet Muenchen,
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
#include "GFDetPlane.h"

#include <assert.h>
#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

ClassImp(GFDetPlane)

GFDetPlane::GFDetPlane(const TVector3& o,
                       const TVector3& u,
                       const TVector3& v,
                       GFAbsFinitePlane* finite)
  :fO(o), fU(u), fV(v), fFinitePlane(finite)
{
  sane();
}


GFDetPlane::GFDetPlane(GFAbsFinitePlane* finite) 
  :fFinitePlane(finite)
{
  // default constructor
  fO.SetXYZ(0.,0.,0.);
  fU.SetXYZ(1.,0.,0.);
  fV.SetXYZ(0.,1.,0.);
  // sane() not needed here
}

GFDetPlane::GFDetPlane(const TVector3& o,
                       const TVector3& n,
                       GFAbsFinitePlane* finite)
  :fO(o), fFinitePlane(finite)
{
  setNormal(n);
}


GFDetPlane::~GFDetPlane(){
  if(fFinitePlane!=NULL) delete fFinitePlane;
}


GFDetPlane::GFDetPlane(const GFDetPlane& rhs) {
  if(rhs.fFinitePlane != NULL) fFinitePlane = rhs.fFinitePlane->clone();
  else fFinitePlane = NULL;
  fO = rhs.fO;
  fU = rhs.fU;
  fV = rhs.fV;
}


GFDetPlane& GFDetPlane::operator=(const GFDetPlane& rhs) {
  if (this == &rhs)
    return *this;
  if(fFinitePlane != NULL) {
    delete fFinitePlane;
  }
  if(rhs.fFinitePlane != NULL){
    fFinitePlane = rhs.fFinitePlane->clone();
  }
  else{
    fFinitePlane = NULL;
  }
  fO = rhs.fO;
  fU = rhs.fU;
  fV = rhs.fV;
  return *this;
}


void 
GFDetPlane::set(const TVector3& o,
                const TVector3& u,
                const TVector3& v)
{
  fO = o;
  fU = u;
  fV = v;
  sane();
}


void 
GFDetPlane::setO(const TVector3& o)
{
  fO = o;
}

void 
GFDetPlane::setO(double X,double Y,double Z)
{
  fO.SetXYZ(X,Y,Z);
}

void 
GFDetPlane::setU(const TVector3& u)
{
  fU = u;
  sane(); // sets fV perpendicular to fU
}

void 
GFDetPlane::setU(double X,double Y,double Z)
{
  fU.SetXYZ(X,Y,Z);
  sane(); // sets fV perpendicular to fU
}

void 
GFDetPlane::setV(const TVector3& v)
{
  fV = v;
  fU = getNormal().Cross(fV);
  fU *= -1.;
  sane();
}

void 
GFDetPlane::setV(double X,double Y,double Z)
{
  fV.SetXYZ(X,Y,Z);
  fU = getNormal().Cross(fV);
  fU *= -1.;
  sane();
}

void 
GFDetPlane::setUV(const TVector3& u,const TVector3& v)
{
  fU = u;
  fV = v;
  sane();
}

void
GFDetPlane::setON(const TVector3& o,const TVector3& n){
  fO = o;
  setNormal(n);
}


TVector3
GFDetPlane::getNormal() const
{
  return fU.Cross(fV);
}

void
GFDetPlane::setNormal(double X,double Y,double Z){
  setNormal( TVector3(X,Y,Z) );
}

void
GFDetPlane::setNormal(const TVector3& n){
  fU = n.Orthogonal();
  fU.SetMag(1.);
  fV = n.Cross(fU);
  fV.SetMag(1.);
}

void GFDetPlane::setNormal(const double& theta, const double& phi){
  setNormal( TVector3(TMath::Sin(theta)*TMath::Cos(phi),TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta)) );
}


TVector2
GFDetPlane::project(const TVector3& x)const
{
  return TVector2(fU*x, fV*x);
}


TVector2 
GFDetPlane::LabToPlane(const TVector3& x)const
{
  return project(x-fO);
}


TVector3
GFDetPlane::toLab(const TVector2& x)const
{
  TVector3 d(fO);
  d += x.X()*fU;
  d += x.Y()*fV;
  return d;
}


TVector3 
GFDetPlane::dist(const TVector3& x)const
{
  return toLab(LabToPlane(x)) - x;
}


void
GFDetPlane::sane(){
  assert(fU!=fV);

  // ensure unit vectors
  fU.SetMag(1.);
  fV.SetMag(1.);

  // check if already orthogonal
  if (fU.Dot(fV) < 1.E-5) return;

  // ensure orthogonal system
  fV = getNormal().Cross(fU);
}


void
GFDetPlane::Print(const Option_t* option) const
{
  std::cout<<"GFDetPlane: "
	   <<"O("<<fO.X()<<","<<fO.Y()<<","<<fO.Z()<<") "
	   <<"u("<<fU.X()<<","<<fU.Y()<<","<<fU.Z()<<") "
	   <<"v("<<fV.X()<<","<<fV.Y()<<","<<fV.Z()<<") "
	   <<"n("<<getNormal().X()<<","<<getNormal().Y()<<","<<getNormal().Z()<<") "
		   <<std::endl;
  std::cout << fFinitePlane << std::endl;
  if(fFinitePlane!=NULL) fFinitePlane->Print(option);
    
}


/*
  I could write pages of comments about correct equality checking for 
  floating point numbers, but: When two planes are as close as 10E-5 cm
  in all nine numbers that define the plane, this will be enough for all
  practical purposes
 */
bool operator== (const GFDetPlane& lhs, const GFDetPlane& rhs){
  if (&lhs == &rhs)
    return true;
  static const double detplaneEpsilon = 1.E-5;
  if(
     fabs( (lhs.fO.X()-rhs.fO.X()) ) > detplaneEpsilon  ||
     fabs( (lhs.fO.Y()-rhs.fO.Y()) ) > detplaneEpsilon  ||
     fabs( (lhs.fO.Z()-rhs.fO.Z()) ) > detplaneEpsilon
     ) return false;
  else if(
	  fabs( (lhs.fU.X()-rhs.fU.X()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fU.Y()-rhs.fU.Y()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fU.Z()-rhs.fU.Z()) ) > detplaneEpsilon
	  ) return false;
  else if(
	  fabs( (lhs.fV.X()-rhs.fV.X()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fV.Y()-rhs.fV.Y()) ) > detplaneEpsilon  ||
	  fabs( (lhs.fV.Z()-rhs.fV.Z()) ) > detplaneEpsilon
	  ) return false;
  return true;
}

bool operator!= (const GFDetPlane& lhs, const GFDetPlane& rhs){
  return !(lhs==rhs);
}


void GFDetPlane::getGraphics(double mesh, double length, TPolyMarker3D **pl,TPolyLine3D** plLine, TPolyLine3D **u, TPolyLine3D **v, TPolyLine3D**n) const {
  *pl = new TPolyMarker3D(21*21,24);
  (*pl)->SetMarkerSize(0.1);
  (*pl)->SetMarkerColor(kBlue);
  int nI=10;
  int nJ=10;
  *plLine = new TPolyLine3D(5);

  {
    TVector3 linevec;
    int i,j;
    i=-1*nI;j=-1*nJ;
    linevec=(fO+(mesh*i)*fU+(mesh*j)*fV);
    (*plLine)->SetPoint(0,linevec.X(),linevec.Y(),linevec.Z());
    i=-1*nI;j=1*nJ;
    linevec=(fO+(mesh*i)*fU+(mesh*j)*fV);
    (*plLine)->SetPoint(0,linevec.X(),linevec.Y(),linevec.Z());
    i=1*nI;j=-1*nJ;
    linevec=(fO+(mesh*i)*fU+(mesh*j)*fV);
    (*plLine)->SetPoint(2,linevec.X(),linevec.Y(),linevec.Z());
    i=1*nI;j=1*nJ;
    linevec=(fO+(mesh*i)*fU+(mesh*j)*fV);
    (*plLine)->SetPoint(1,linevec.X(),linevec.Y(),linevec.Z());
    i=-1*nI;j=-1*nJ;
    linevec=(fO+(mesh*i)*fU+(mesh*j)*fV);
    (*plLine)->SetPoint(4,linevec.X(),linevec.Y(),linevec.Z());
  }
  for (int i=-1*nI;i<=nI;++i){
    for (int j=-1*nJ;j<=nJ;++j){
      TVector3 vec(fO+(mesh*i)*fU+(mesh*j)*fV);
      int id=(i+10)*21+j+10;
      (*pl)->SetPoint(id,vec.X(),vec.Y(),vec.Z());
    }
  }


  *u = new TPolyLine3D(2);
  (*u)->SetPoint(0,fO.X(),fO.Y(),fO.Z());
  (*u)->SetPoint(1,fO.X()+length*fU.X(),fO.Y()+length*fU.Y(),fO.Z()+length*fU.Z());
  (*u)->SetLineWidth(2);
  (*u)->SetLineColor(kGreen);


  *v = new TPolyLine3D(2);
  (*v)->SetPoint(0,fO.X(),fO.Y(),fO.Z());
  (*v)->SetPoint(1,fO.X()+length*fV.X(),fO.Y()+length*fV.Y(),fO.Z()+length*fV.Z());
  (*v)->SetLineWidth(2);
  (*v)->SetLineColor(kRed);

  if(n!=NULL){
    *n = new TPolyLine3D(2);
    TVector3 _n=getNormal();
    (*n)->SetPoint(0,fO.X(),fO.Y(),fO.Z());
    (*n)->SetPoint(1,fO.X()+length*_n.X(),fO.Y()+length*_n.Y(),fO.Z()+length*_n.Z());
    (*n)->SetLineWidth(2);
    (*n)->SetLineColor(kBlue);
  }
}


double GFDetPlane::distance(const TVector3& point) const {
  // |(point - fO)*(fU x fV)|
  return fabs( (point.X()-fO.X()) * (fU.Y()*fV.Z() - fU.Z()*fV.Y()) +
               (point.Y()-fO.Y()) * (fU.Z()*fV.X() - fU.X()*fV.Z()) +
               (point.Z()-fO.Z()) * (fU.X()*fV.Y() - fU.Y()*fV.X()));
}

double GFDetPlane::distance(double x, double y, double z) const {
  // |(point - fO)*(fU x fV)|
  return fabs( (x-fO.X()) * (fU.Y()*fV.Z() - fU.Z()*fV.Y()) +
               (y-fO.Y()) * (fU.Z()*fV.X() - fU.X()*fV.Z()) +
               (z-fO.Z()) * (fU.X()*fV.Y() - fU.Y()*fV.X()));
}


TVector2 GFDetPlane::straightLineToPlane (const TVector3& point, const TVector3& dir) const {
  TVector3 dirNorm(dir.Unit());
  TVector3 normal = getNormal();
  double dirTimesN = dirNorm*normal;
  if(fabs(dirTimesN)<1.E-6){//straight line is parallel to plane, so return infinity
    return TVector2(1.E100,1.E100);
  }
  double t = 1./dirTimesN * ((fO-point)*normal);
  return project(point - fO + t * dirNorm);
}


void GFDetPlane::reset() {
  fO.SetXYZ(0.,0.,0.);
  fU.SetXYZ(1.,0.,0.);
  fV.SetXYZ(0.,1.,0.);
  if(fFinitePlane!=NULL) {
    delete fFinitePlane;
    fFinitePlane = NULL;
  }
}

