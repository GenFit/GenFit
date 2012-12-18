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
#ifndef RKTRACKREP_H
#define RKTRACKREP_H

/*
  Track parametrization adapted from the COMPASS experimt's analysis 
  software Phast. The state vector is (X,Y,dX/dZ,dY/dZ,p/|vec{p}|)
 */

#include "GFAbsTrackRep.h"

class RKTrackRepXY : public GFAbsTrackRep {
public:

  // Constructors/Destructors ---------
  RKTrackRepXY();
  RKTrackRepXY(const TVector3& pos,
	     const TVector3& mom,
	     const TVector3& poserr,
	     const TVector3& momerr,
	     //const double& q,
	     const int& PDGCode);

  virtual ~RKTrackRepXY();


  virtual GFAbsTrackRep* clone() const {return new RKTrackRepXY(*this);}
  virtual GFAbsTrackRep* prototype()const{return new RKTrackRepXY();}

  virtual double extrapolate(const GFDetPlane&, 
			   TMatrixT<double>& statePred,
			   TMatrixT<double>& covPred);

  
  virtual double extrapolateToPoint(const TVector3& pos,
				 TVector3& poca,
				 TVector3& dirInPoca);
  
  virtual double extrapolateToLine(const TVector3& point1,
				 const TVector3& point2,
				 TVector3& poca,
				 TVector3& dirInPoca,
				 TVector3& poca_onwire);
  

  virtual TVector3 getPos(const GFDetPlane&);
  virtual TVector3 getMom(const GFDetPlane&);
  virtual void getPosMom(const GFDetPlane&,TVector3& pos,TVector3& mom);

  virtual double getCharge()const {return fCharge;}

  void switchDirection(){fDirection = (!fDirection);}

  void setPDG(int);

private:
  bool fDirection;
  int fPdg;
  double fMass;
  double fCharge;

  TMatrixT<double> cov15to25(const TMatrixT<double>& cov15) const;
  TMatrixT<double> cov25to15(const TMatrixT<double>& cov25) const;

  double getStep(const double& zFinal, const double& distance, const TVector3& pos, const TVector3& mom,double& XX0, double& dP);
  void addNoise(double len,const TMatrixT<double>& s,TMatrixT<double>& cov15);
  bool RKutta (double* SU,double* VO, double& Path) const;
  double Extrap( double Z, const double& zFrom,double& zOut, const TMatrixT<double>& stateIn,TMatrixT<double>& stateOut, const TMatrixT<double>& covIn,TMatrixT<double>& covOut) const;
 public:
  ClassDef(RKTrackRepXY,2)


};


#endif

