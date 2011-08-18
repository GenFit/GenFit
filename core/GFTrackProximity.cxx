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
#include "TVector3.h"
#include <iostream>
#include "GFTrackProximity.h"
#include "GFTrack.h"
#include "GFAbsTrackRep.h"
using namespace std;

double trkDist(GFAbsTrackRep* rep1, GFAbsTrackRep* rep2){

  TVector3 d=rep1->getPos()-rep2->getPos();
  //cout<<"trkDist"<<d.Mag()<<endl;
  return d.Mag();
}

TVector3 trackProximity(GFTrack *trk1, GFTrack *trk2){

  GFAbsTrackRep* rep1=trk1->getCardinalRep();
  GFAbsTrackRep* rep2=trk2->getCardinalRep();
  return trackProximity(rep1,rep2);
}

TVector3 trackProximity(GFAbsTrackRep* rep1, GFAbsTrackRep* rep2){
  // TODO: make accuracy configurable!
  // a simple newtonian search.
  double h=0.01;
  double m1=999999;
  double s1=0;
  int steps=0;
  double oldd=999999;
  double d=99999;
  //  TVector3 pos,mom,poserr,momerr;
 
  
  //cout<<"GFTrackProximity start ########################################################################################################################"<<endl;


  
  while(oldd>d && steps<100){
    
    //cout<<"GFTrackProximity::steps:"<<steps<<endl;
    TVector3 point1,dir1;    
    TVector3 point2, norm2;
    GFDetPlane plane1,plane2;
    rep1->stepalong(s1,point1,dir1);
    plane1.setO(point1);
    //cout<<"track1 ";point1.Print();
    plane1.setNormal(dir1);
    rep1->GFAbsTrackRep::extrapolate(plane1);
    rep2->extrapolateToPoint(point1,point2,norm2);
    plane2.setO(point1);
    plane2.setNormal(norm2);
    rep2->GFAbsTrackRep::extrapolate(plane2);

    oldd=d;
    ++steps;
    d=trkDist(rep1,rep2);
    //cout<<"dist"<<d<<endl;
    //----------------------------------------------------------------------
    rep1->stepalong(h,point1,dir1);
    //cout<<"track1 ";point1.Print();
    plane1.setO(point1);
    plane1.setNormal(dir1);
    rep1->GFAbsTrackRep::extrapolate(plane1);
    rep2->extrapolateToPoint(point1,point2,norm2);
    plane2.setO(point1);
    plane2.setNormal(norm2);
    rep2->GFAbsTrackRep::extrapolate(plane2);

    double d1=trkDist(rep1,rep2);
    //cout<<"dist"<<d1<<endl;
    m1=(d1-d)/h;
    if(TMath::Abs(m1)<1E-4)break;
    s1=-d/m1;
    // limit stepsize to max 180cm
    if(TMath::Abs(s1)>200)s1= s1>0 ? 180 : -180;
    //std:://cout<<"d="<<d<<"  d1="<<d1<<"  m1="<<m1<<"  s1="<<s1<<std::endl;
    
  }
  TVector3 point1,dir1;
  GFDetPlane plane1, plane2;
  TVector3 point2, norm2;
  rep1->stepalong(s1,point1,dir1);
  //cout<<"track1 ";point1.Print();
  plane1.setO(point1);
  plane1.setNormal(dir1);
  rep1->GFAbsTrackRep::extrapolate(plane1);
  rep2->extrapolateToPoint(point1,point2,norm2);
  plane2.setO(point1);
  plane2.setNormal(norm2);
  rep2->GFAbsTrackRep::extrapolate(plane2);
  
  // three points -> fit a parabola -> minimum
  h=1.;
  double d0=trkDist(rep1,rep2);
  //cout<<"dist"<<d0<<endl;
  rep1->stepalong(h,point1,dir1);
  //cout<<"track1 ";point1.Print();
  plane1.setO(point1);
  plane1.setNormal(dir1);
  rep1->GFAbsTrackRep::extrapolate(plane1);
  rep2->extrapolateToPoint(point1,point2,norm2);
  plane2.setO(point1);
  plane2.setNormal(norm2);
  rep2->GFAbsTrackRep::extrapolate(plane2);
  double d1=trkDist(rep1,rep2);
  //cout<<"dist"<<d1<<endl;
  rep1->stepalong(-2.*h,point1,dir1);
  //cout<<"track1 ";point1.Print();
  plane1.setO(point1);
  plane1.setNormal(dir1);
  rep1->GFAbsTrackRep::extrapolate(plane1);
  rep2->extrapolateToPoint(point1,point2,norm2);
  plane2.setO(point1);
  plane2.setNormal(norm2);
  rep2->GFAbsTrackRep::extrapolate(plane2);
  double d2=trkDist(rep1,rep2);
  //cout<<"dist"<<d2<<endl;
  
  double s=(d2-d1)*h;
  double denom=2.*(d1+d2-2*d0);
  if(denom==0){
    TVector3 posA=rep1->GFAbsTrackRep::getPos();
    TVector3 posB=rep2->GFAbsTrackRep::getPos();
    TVector3 poca=(posB+posA)*0.5;
    return poca;
    //    return trkDist(rep1,rep2);
  }
  s/=denom;

  //std:://cout<<"d0="<<d0<<"  d1="<<d1<<"  d2="<<d2<<"   s="<<s<<std::endl;

  rep1->stepalong(s+h,point1,dir1);
  //cout<<"track1 ";point1.Print();
  plane1.setO(point1);
  plane1.setNormal(dir1);
  rep1->GFAbsTrackRep::extrapolate(plane1);
  rep2->extrapolateToPoint(point1,point2,norm2);
  plane2.setO(point1);
  plane2.setNormal(norm2);
  rep2->GFAbsTrackRep::extrapolate(plane2);
  //double d3=trkDist(rep1,rep2);
  //cout<<"dist"<<d3<<endl;
  //cout<<"GFTrackProximity end ########################################################################################################################"<<endl;
  TVector3 posA=rep1->GFAbsTrackRep::getPos();
  TVector3 posB=rep2->GFAbsTrackRep::getPos();

  TVector3 poca=(posB+posA)*0.5;
  return poca;
  //return trkDist(rep1,rep2);
}


