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
//      4th order Runge Kutta algorithm after Nystrom
//      integrates a second order ordinary differentiel equation of type
//      u''=g(u,u'); u=(u1,u2,...un); u\element{R}
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef NYSTROM_HH
#define NYSTROM_HH

// Base Class Headers ----------------


// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op
#include "TVectorT.h"

// Collaborating Class Declarations --
class AbsNystromEQM;


class Nystrom {
public:

  // Constructors/Destructors ---------
  Nystrom(AbsNystromEQM* eqm);
  ~Nystrom(){}
  
  // Modifiers -----------------------
  void setFunction(AbsNystromEQM* eqm){_eqm=eqm;}
  void setAccuracy(double a){_acc=a;}
  void setAdaptive(bool flag=true){_adaptive=flag;}
  
  // Operations ----------------------
  void step(const TVectorT<double>& u, 
	    const TVectorT<double>& uprime,
	    const TVectorT<double>& par,
	    TVectorT<double>& newu, TVectorT<double>& newuprime,
	    double h);

  double adaptiveStep(const TVectorT<double>& u, 
		      const TVectorT<double>& uprime, 
		      const TVectorT<double>& par,
		      TVectorT<double>& newu,
		      TVectorT<double>& newuprime,
		      double& stepdone,
		      double h,double sign, double reststep);  // return new h


  // returns length of propagation (in 3D!)
  double propagate(double start, double end,
		 const TVectorT<double>& u, 
		 const TVectorT<double>& uprime,
		 const TVectorT<double>& par,
		 TVectorT<double>& newu, TVectorT<double>& newuprime);
private:
  
  // Private Data Members ------------

  //  Equation of Motion:
  AbsNystromEQM* _eqm;
  double _acc; // desired accuracy
  bool _adaptive; // use adaptive stepping? default=true;

  // function pointer to system function u''=g(u,u')
  // TVectorT<double>(*_g)(const TVectorT<double>&u,const TVectorT<double>& uprim, const TVectorT<double>& par);
  
  // Private Methods -----------------

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
