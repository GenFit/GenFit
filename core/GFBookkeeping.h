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



#ifndef GFBOOKKEEPING_H
#define GFBOOKKEEPING_H

#include"TObject.h"
#include"TMatrixT.h"
#include<vector>
#include<cassert>
#include<iostream>
#include<utility>
#include<map>
#include"GFDetPlane.h"

class GFBookkeeping : public TObject {
 private:

  //the string keys will in general be different, so this cant
  //be unified to one container
  std::map<std::string, TMatrixT<double>* > fMatrices;
  std::map<std::string, GFDetPlane* > fPlanes;
  /* this is a work-around: we want to save doubles, but ROOT has problems
   * with TObjects that contain map<string,double*>. We take a 1x1 matrix
   * as a work-around to hold the double internally */
  std::map<std::string, TMatrixT<double>* > fNumbers;
  std::vector< unsigned int > fFailedHits;
  int fNhits;

 public:
  void reset();
  void setNhits(int n){fNhits=n; reset();}

  void bookMatrices(std::string key);
  void bookGFDetPlanes(std::string key);
  void bookNumbers(std::string key,double val=0.);

  void setMatrix(std::string key,unsigned int index,const TMatrixT<double>& mat);
  void setDetPlane(std::string key,unsigned int index,const GFDetPlane& pl);
  void setNumber(std::string key,unsigned int index, const double& num);

  bool getMatrix(std::string key, unsigned int index, TMatrixT<double>& mat) const;
  bool getDetPlane(std::string key, unsigned int index, GFDetPlane& pl)  const;
  bool getNumber(std::string key, unsigned int index, double& num) const;

  std::vector< std::string > getMatrixKeys() const;
  std::vector< std::string > getGFDetPlaneKeys() const;
  std::vector< std::string > getNumberKeys() const;

  void addFailedHit(unsigned int);
  unsigned int hitFailed(unsigned int);
  unsigned int getNumFailed();

  GFBookkeeping(){fNhits=-1;}
  GFBookkeeping(const GFBookkeeping&);
  virtual ~GFBookkeeping(){clearAll();}

  void clearAll();
  void clearFailedHits();

  void Print(const Option_t* = "") const;

 private:
  //protect from call of not yet defined assignement operator
  GFBookkeeping& operator=(const GFBookkeeping&){
    return *this;
  }

 public:
  ClassDef(GFBookkeeping,2)

};

#endif
