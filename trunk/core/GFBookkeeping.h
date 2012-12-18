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

#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <vector>
#include <iostream>
#include <map>
#include "GFDetPlane.h"


enum GFBKKey {
  // forward info ------
  GFBKKey_fSt,
  GFBKKey_fCov,
  GFBKKey_fUpSt,
  GFBKKey_fUpCov,
  GFBKKey_fPl,
  GFBKKey_fAuxInfo,
  GFBKKey_fExtLen,
  // add further forward keys here

  // backward info ------
  GFBKKey_bSt=1000,
  GFBKKey_bCov,
  GFBKKey_bUpSt,
  GFBKKey_bUpCov,
  GFBKKey_bPl,
  GFBKKey_bAuxInfo,
  GFBKKey_bExtLen,
  // add further backward keys here

  // other info ------
  GFBKKey_dafWeight=2000
  // add further other keys here
};


class GFBookkeeping : public TObject {
 private:

  //the string keys will in general be different, so this can't
  //be unified to one container
  std::map<GFBKKey, std::vector<TVectorD> > fVectors;
  std::map<GFBKKey, std::vector<TMatrixD> > fMatrices;
  std::map<GFBKKey, std::vector<TMatrixDSym> > fSymMatrices;
  std::map<GFBKKey, std::vector<GFDetPlane> > fPlanes;
  std::map<GFBKKey, std::vector<double> > fNumbers;
  std::vector< unsigned int > fFailedHits;
  int fNhits;

 public:
  /** @brief Keep all keys and clear only the data.
   * Resize the data vectors to fNhits.
   * Also clears fFailedHits.
   */
  void reset();

  /** @brief Set number of hits, clear data but keep keys.
   */
  void setNhits(int n){fNhits=n; reset();}

  void bookVectors(const GFBKKey& key);
  void bookMatrices(const GFBKKey& key);
  void bookSymMatrices(const GFBKKey& key);
  void bookGFDetPlanes(const GFBKKey& key);
  void bookNumbers(const GFBKKey& key,double val=0.);

  void setVector(const GFBKKey& key,unsigned int index,const TVectorD& mat);
  void setMatrix(const GFBKKey& key,unsigned int index,const TMatrixD& mat);
  void setSymMatrix(const GFBKKey& key,unsigned int index,const TMatrixDSym& mat);
  void setDetPlane(const GFBKKey& key,unsigned int index,const GFDetPlane& pl);
  void setNumber(const GFBKKey& key,unsigned int index, const double& num);

  const TVectorD&    getVector(const GFBKKey& key, unsigned int index) const;
  const TMatrixD&    getMatrix(const GFBKKey& key, unsigned int index) const;
  const TMatrixDSym& getSymMatrix(const GFBKKey& key, unsigned int index) const;
  const GFDetPlane&  getDetPlane(const GFBKKey& key, unsigned int index)  const;
  double             getNumber(const GFBKKey& key, unsigned int index) const;

  std::vector<GFBKKey> getVectorKeys() const;
  std::vector<GFBKKey> getMatrixKeys() const;
  std::vector<GFBKKey> getSymMatrixKeys() const;
  std::vector<GFBKKey> getGFDetPlaneKeys() const;
  std::vector<GFBKKey> getNumberKeys() const;

  void addFailedHit(unsigned int);
  /** @brief How often did hit nr. iHit fail
   */
  unsigned int hitFailed(unsigned int iHit);
  unsigned int getNumFailed();

  GFBookkeeping(){fNhits=-1;}
  GFBookkeeping(const GFBookkeeping&);

  /** @brief clear fVectors, fMatrices, fSymMatrices, fPlanes, fNumbers
   */
  void clearAll();

  /** @brief clear fFailedHits
   */
  void clearFailedHits();

  void Print(const Option_t* = "") const;

 private:
  //protect from call of not yet defined assignment operator
  GFBookkeeping& operator=(const GFBookkeeping&){
    return *this;
  }

 public:
  ClassDef(GFBookkeeping,4)

};


#endif
