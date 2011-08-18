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
//      Implementation of class LSLEQM
//      see LSLEQM.hh for details
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
#include "LSLEQM.h"

// C/C++ Headers ----------------------
#include <iostream>
#include "TMath.h"
// Collaborating Class Headers --------
#include "GFAbsBField.h"

// Class Member definitions -----------


TVectorT<double> 
LSLEQM::eval(const TVectorT<double>&u,
	     const TVectorT<double>& uprim,
	     const TVectorT<double>& par) {
  // default values in case no bfield is given
  double Bx=0.; double By=0.;
  double Bz=2.;
  if(u[2]>200.)Bz=0.;
  if(u[2]>350. && u[2]<450.)By=0.5;

  TVector3 uu(u[0],u[1],u[2]);
  if(_field!=NULL){ 
    //std::cout<<"using fieldmap"<<std::endl;
    TVector3 B=_field->get(uu);
    //B field is now in kGauss, so convert to Tesla
    B*=0.1;
    Bx=B.X();By=B.Y();Bz=B.Z();;
  }
  
  //if(uprim[0]!=1.)
  //  std::cout<<"uprim("<<uprim[0]<<","<<uprim[1]<<")"<<std::endl;

// limit possible slopes
  double xprim=uprim[0];
  double yprim=uprim[1];
  if(TMath::Abs(xprim)>1E4)xprim<0 ? xprim=-1E4 : xprim=1E4;
  if(TMath::Abs(yprim)>1E4)yprim<0 ? yprim=-1E4 : yprim=1E4;
  double xprim2=xprim*xprim;
  double yprim2=yprim*yprim;
  double dsdz=TMath::Sqrt(1.+xprim2+yprim2);
  double kappP=0.00299792458*par[0]; // [GeV/c T^-1 cm^-1]
  TVectorT<double> result(u);

  result[0]=kappP*dsdz*(xprim*yprim*Bx-(1+xprim2)*By+yprim*Bz);
  result[1]=kappP*dsdz*(-xprim*yprim*By+(1+yprim2)*Bx-xprim*Bz);
  result[2]=0;
  return result;
}
