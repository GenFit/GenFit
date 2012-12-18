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

#include "GFAbsTrackRep.h"

GFAbsTrackRep::GFAbsTrackRep() : fDimension(5),fState(5), fCov(5), fChiSqu(0), fForwardChiSqu(0), fNdf(0), fStatusFlag(0), fFirstState(5), fFirstCov(5), fLastState(5), fLastCov(5), fXX0(-1.)
{
}

GFAbsTrackRep::GFAbsTrackRep(int dim) : fDimension(dim), fState(dim), fCov(dim), fChiSqu(0), fForwardChiSqu(0), fNdf(0), fStatusFlag(0), fFirstState(dim), fFirstCov(dim), fLastState(dim), fLastCov(dim), fXX0(-1.)
{
}

GFAbsTrackRep::~GFAbsTrackRep() {}

double GFAbsTrackRep::extrapolate(const GFDetPlane& plane){
  TVectorD statePred(fDimension);
  TMatrixDSym covPred(fDimension);
  double retVal = extrapolate(plane,statePred,covPred);
  setData(statePred,plane,&covPred);
  return retVal;
}

//default implentation might be overwritten, please see the doxy docu
double GFAbsTrackRep::extrapolate(const GFDetPlane& plane, TVectorD& statePred){
  TMatrixDSym cov(fDimension);
  return extrapolate(plane,statePred,cov);
}

void GFAbsTrackRep::Abort(std::string method){
  std::cerr << method <<  " as implemented in " << __FILE__ 
	    << " was called. This means that this feature was used "
	    << "in a track rep which didnt overwrite this method. "
	    << std::endl << "C++ throw;" << std::endl;
  //system call abort
  throw;
}

double GFAbsTrackRep::extrapolateToPoint(const TVector3& point,
				    TVector3& poca,
				    TVector3& normVec){
  Abort("extrapolateToPoint()");
  return 0.;
}

double GFAbsTrackRep::extrapolateToLine(const TVector3& point1,
									const TVector3& point2,
									TVector3& poca,
									TVector3& normVec,
									TVector3& poca_onwire){
  Abort("extrapolateToLine()");
  return 0.;
}
  

double GFAbsTrackRep::stepalong(double h,
                  TVector3& point,
                  TVector3& dir){
  Abort("stepalong()");
  return 0.;
}


int GFAbsTrackRep::getPDG(){
  Abort("getPDG()");
  return 0;
}


void GFAbsTrackRep::getPosMomCov(const GFDetPlane& pl, TVector3& pos, TVector3& mom, TMatrixDSym& cov){
  Abort("getPosMomCov()");
}

void GFAbsTrackRep::setPosMomCov(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov){
  Abort("setPosMomCov()");
}

void
GFAbsTrackRep::reset(){
  std::cout<<"GFAbsTrackRep::reset"<<std::endl;
  fRefPlane.reset();
  fState.Zero();
  fCov.Zero();
  fFirstState.Zero();
  fFirstCov.Zero();
  fLastState.Zero();
  fLastCov.Zero();
}

void
GFAbsTrackRep::Print(const Option_t* option) const {
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout<<"GFAbsTrackRep::Parameters at reference plane ";
  fRefPlane.Print(option);
  std::cout<<"GFAbsTrackRep::State"<<std::endl;
  fState.Print(option);
  std::cout<<"GFAbsTrackRep::Covariances"<<std::endl;
  fCov.Print(option);
  std::cout<<"GFAbsTrackRep::chi^2"<<std::endl;
  std::cout<<fChiSqu<<std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
}



ClassImp(GFAbsTrackRep)
