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
//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class Nystrom
//      see Nystrom.hh for details
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
#include "Nystrom.h"
#include <cmath>
// C/C++ Headers ----------------------
#include <iostream>
#include "assert.h"
#include "TMath.h"

// Collaborating Class Headers --------
#include "AbsNystromEQM.h"

// Class Member definitions -----------


Nystrom::Nystrom(AbsNystromEQM* eqm)
  : _eqm(eqm), _acc(1E-2), _adaptive(true)
{}

void
Nystrom::step(const TVectorT<double>& u, const TVectorT<double>& uprime, 
	      const TVectorT<double>& par,
	      TVectorT<double>& newu, TVectorT<double>& newuprime,
	      double h)
{
  assert(h!=0);
  double sign= h > 0 ? 1. : -1.;
  double habs=sign*h;
  double hh=habs*0.5;
  double h2=h*h;
  TVectorT<double> k1=_eqm->eval(u,uprime,par);
  TVectorT<double> k2=_eqm->eval(u+hh*uprime+h2/8.*k1,uprime+hh*k1,par);
  TVectorT<double> k3=_eqm->eval(u+hh*uprime+h2/8.*k2,uprime+hh*k2,par);
  TVectorT<double> k4=_eqm->eval(u+habs*uprime+h2/2.*k3,uprime+habs*k3,par); 
  
  TVectorT<double> ksum=k1+k2+k3;
  newu=u+h*uprime+h2/6.*ksum;
  newuprime=uprime+(h/6.)*(ksum+k4+k2+k3);
  assert(TMath::Abs(newuprime[2]-1.)<1E-4);
}


// this will return the new h!
double
Nystrom::adaptiveStep(const TVectorT<double>& u, 
		      const TVectorT<double>& uprime, 
		      const TVectorT<double>& par,
		      TVectorT<double>& newu,
		      TVectorT<double>& newuprime,
		      double& stepdone,
		      double h, double sign, double reststep)
{
  
  //std::cout<<"Starting adaptive step:"<<std::endl;
  double d0=_acc*h; // desired accuracy
  double d=1.1*d0;
  double newh=h;
  double _maxstep=0.1;
  

  assert(reststep>=0);

  TVectorT<double> newu1(u);
  TVectorT<double> newuprime1(uprime);
  TVectorT<double> newu2(u);
  TVectorT<double> newuprime2(uprime);
  TVectorT<double> utmp(u);
  TVectorT<double> uprimetmp(uprime);
  TVectorT<double> delta(u);
  //newh=0.01;
  int n=delta.GetNrows();
  while(d>d0){

    if(newh>=reststep)newh=reststep;
    
    //std::cout<<"newh="<<newh<<"  d="<<d<<std::endl;
    // make one step
    step(u,uprime,par,
	 newu1,newuprime1,newh*sign);

  
    // make two half-steps
    step(u,uprime,par,
	 utmp,uprimetmp,newh*0.5*sign);
    step(utmp,uprimetmp,par,
	 newu2,newuprime2,newh*0.5*sign);
  
    stepdone=newh*sign;

    // calculate delta:
    delta=newu1-newu2;
    //delta.Print();

    // get biggest contribution:
    int index=0;
    d=0;
    for(int i=0; i<n; ++i){
      if(TMath::Abs(delta[i])>d){
	d=TMath::Abs(delta[i]);
	index=i;
      }
    }
    d0=TMath::Abs(_acc*h*uprime[index]);

    //std::cout<<"d="<<d<<"   d0="<<d0
    //     <<"   zdif="<<newu1[2]-newu2[2]
    //	     <<"   newh="<<newh<<std::endl;

    assert(d>=0);
    if(d>0){
      double frac=TMath::Abs(d0/d);
      newh=0.9*TMath::Abs(newh*(TMath::Power(frac,0.25)));
      if(newh<0.00001){
	newh=0.00001;
	break;
      }
      if(newh>_maxstep)newh=_maxstep;
    }
    else {
      //std::cout<<"small dev! newh="<<newh<<std::endl;
      newh=2*newh;
      if(newh>_maxstep)newh=_maxstep;
      break;
    }
    //break;
  }
  
  //std::cout<<"completed step d="<<d
  // 	   <<"  d0="<<d0
  //	   <<"  newh="<<newh
  //	   <<"  stepdone="<<stepdone<<std::endl;
  newu=newu1;//+0.06666666667*delta;
  newuprime=newuprime1;
  assert(newh>0);
  return newh;
}


// returns length
double 
Nystrom::propagate(double start, double end,
		   const TVectorT<double>& u, 
		   const TVectorT<double>& uprime,
		   const TVectorT<double>& par,
		   TVectorT<double>& newu, TVectorT<double>& newuprime)
{
  if(TMath::Abs(start-end)<1E-8){
    newu=u;
    newuprime=uprime;
    return 0;
  }

  double length=0;

  double sign= start < end ? 1. : -1.;
  double h=_acc; // TODO: This should be more dynamic!
  double s=start;
  TVectorT<double> utmp(u);
  TVectorT<double> uprimetmp(uprime);
  unsigned int steps=0;
  unsigned int maxsteps=1000000;
  double z0=u[2];
  double z=z0;
  assert(fabs(z0-start)<1E-8);

  while(TMath::Abs(s-end)>h && steps<maxsteps){

    ++steps;
    
    double reststep=TMath::Abs(s-end);
    double newh=h;
    double stepdone=sign*h;
    if(_adaptive){
      newh=adaptiveStep(utmp,uprimetmp,par,
			newu,newuprime,stepdone,h,sign,reststep);
    }
    else {
      step(utmp,uprimetmp,par,
	   newu,newuprime,sign*h);
    }
    //if(fabs(newh-h)>1E-4)std::cout<<"h="<<h<<"  newh="<<newh<<std::endl;
    double dx=h*uprimetmp[0];
    double dy=h*uprimetmp[1];

    length+=sqrt(dx*dx+dy*dy+h*h);


    s+=stepdone;
    z+=stepdone;

    h=newh;
    utmp=newu;
    uprimetmp=newuprime;

    // todo: remove hard limits
    if(fabs(utmp[0])>1000 || fabs(utmp[1])>1000){
      std::cerr<<"Nystrom:: outside of range (-1000,1000)!"<<std::endl;
      return length;
    }
    
  }

 
  

  if(steps>=maxsteps)std::cerr<<"Nystrom:: too many steps (>"<<maxsteps<<") aborting."<<std::endl;
  // do the rest step
  h=end-s;
  if(h>1E-8){
    step(utmp,uprimetmp,par,
	 newu,newuprime,h);
  }
  double dx=h*uprimetmp[0];
  double dy=h*uprimetmp[1];
  
  length+=sqrt(dx*dx+dy*dy+h*h);
  return length;
  //double avstep=fabs(z0-z)/(double)steps;
  //std::cout<<"z0="<<z0
  //	   <<"   z="<<z
  //   <<"   newu[2]="<<newu[2]
  //   <<"   avstep="<<avstep
  //   <<"   signe="<<sign<<std::endl;

}
