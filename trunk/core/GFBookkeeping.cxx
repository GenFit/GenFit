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

#include "GFBookkeeping.h"
#include "GFException.h"


GFBookkeeping::GFBookkeeping(const GFBookkeeping& bk) {
  fNhits = bk.fNhits;
  fVectors = bk.fVectors;
  fMatrices = bk.fMatrices;
  fSymMatrices = bk.fSymMatrices;
  fNumbers = bk.fNumbers; 
  fPlanes = bk.fPlanes; 
  fFailedHits = bk.fFailedHits;
}

// Use a custom streamer, the auto-generated one is prohibitively slower (root 5.34).
void GFBookkeeping::Streamer(TBuffer &R__b)
{
   // Stream an object of class GFBookkeeping.
   if (R__b.IsReading()) {

     //     Version_t R__v = R__b.ReadVersion(); 
     TObject::Streamer(R__b); 

     clearAll();

     R__b >> fNhits;

     unsigned int keyInt;
     GFBKKey key;
     unsigned int nkeys;
     TVectorD vec;
     TMatrixD mat;
     TMatrixDSym symmat;
     GFDetPlane pl;

     {//reading vectors
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         R__b >> keyInt;
         key = (GFBKKey)keyInt;
         bookVectors(key);
         for(int j=0;j<fNhits;++j){
           vec.Streamer(R__b);
           setVector(key,j,vec);
         }
       }
     }//done reading vectors
     {//reading matrices
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         R__b >> keyInt;
         key = (GFBKKey)keyInt;
         bookMatrices(key);
         for(int j=0;j<fNhits;++j){
           mat.Streamer(R__b);
           setMatrix(key,j,mat);
         }
       }
     }//done reading matrices
     {//reading symmetric matrices
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         R__b >> keyInt;
         key = (GFBKKey)keyInt;
         bookSymMatrices(key);
         for(int j=0;j<fNhits;++j){
           symmat.Streamer(R__b);
           setSymMatrix(key,j,symmat);
         }
       }
     }//done reading matrices
     {//reading planes
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         R__b >> keyInt;
         key = (GFBKKey)keyInt;
         bookGFDetPlanes(key);
         for(int j=0;j<fNhits;++j){
           pl.Streamer(R__b);
           setDetPlane(key,j,pl);
         }
       }
     }//done reading planes
     {//reading numbers
       R__b >> nkeys;
       for(unsigned int i=0;i<nkeys;++i){
         R__b >> keyInt;
         key = (GFBKKey)keyInt;
         bookNumbers(key);
         for(int j=0;j<fNhits;++j){
           double d;
           R__b >> d;
           setNumber(key,j,d);
         }
       }
     }//done reading numbers
     {//read failed hits
       clearFailedHits();
       unsigned int nFailedHits;
       R__b >> nFailedHits;
       fFailedHits.reserve(nFailedHits);
       for(unsigned int i=0;i<nFailedHits;++i){
         unsigned int aFailedHit;
         R__b >> aFailedHit;
         fFailedHits[i] = aFailedHit;
       }
     }//done reading failed hits
   } else {
     //     R__b.WriteVersion(GFBookkeeping::IsA()); 
     TObject::Streamer(R__b); 

     //write number of hits
     R__b << fNhits;

     std::vector<GFBKKey> keys;
     {//save vectors
       keys = getVectorKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         R__b << (unsigned int)(keys[i]);
         for(int j=0;j<fNhits;++j){
           ((fVectors[keys[i]])[j]).Streamer(R__b);
         }
       }
     }
     {//save matrices
       keys = getMatrixKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         R__b << (unsigned int)(keys[i]);
         for(int j=0;j<fNhits;++j){
           ((fMatrices[keys[i]])[j]).Streamer(R__b);
         }
       }
     }
     {//save symmetric matrices
       keys = getSymMatrixKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         R__b << (unsigned int)(keys[i]);
         for(int j=0;j<fNhits;++j){
           ((fSymMatrices[keys[i]])[j]).Streamer(R__b);
         }
       }
     }
     keys.clear();
     {//save GFDetPlanes
       keys = getGFDetPlaneKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         R__b << (unsigned int)(keys[i]);
         for(int j=0;j<fNhits;++j){
           ((fPlanes[keys[i]])[j]).Streamer(R__b);
         }
       }
     }//done saving GFDetPlanes
     keys.clear();
     {//save numbers    
       keys = getNumberKeys();
       R__b << (unsigned int)(keys.size());
       for(unsigned int i=0;i<keys.size();++i){
         R__b << (unsigned int)(keys[i]);
         for(int j=0;j<fNhits;++j){
           R__b << (fNumbers[keys[i]])[j];
         }
       }
     }//done saving numbers 
     {//save failedHits
       R__b << ((unsigned int) fFailedHits.size());
       for(unsigned int i=0;i<fFailedHits.size();++i){
         R__b << fFailedHits[i];
       }
     }//done saving failed Hits    
   }
}


void GFBookkeeping::bookVectors(const GFBKKey& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<GFBKKey, std::vector<TVectorD> >::const_iterator it;
  it = fVectors.find(key);
  if(it != fVectors.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookVectors()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fVectors[key].resize(fNhits); // valgrind detects a memory leak here, but it is probably only due to stl optimizations: http://stackoverflow.com/questions/1379564/c-tiny-memory-leak-with-stdmap
}

void GFBookkeeping::bookMatrices(const GFBKKey& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<GFBKKey, std::vector<TMatrixD> >::const_iterator it;
  it = fMatrices.find(key);
  if(it != fMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookMatrices()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fMatrices[key].resize(fNhits); // valgrind detects a memory leak here, but it is probably only due to stl optimizations: http://stackoverflow.com/questions/1379564/c-tiny-memory-leak-with-stdmap
}

void GFBookkeeping::bookSymMatrices(const GFBKKey& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<GFBKKey, std::vector<TMatrixDSym> >::const_iterator it;
  it = fSymMatrices.find(key);
  if(it != fSymMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookSymMatrices()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fSymMatrices[key].resize(fNhits); // valgrind detects a memory leak here, but it is probably only due to stl optimizations: http://stackoverflow.com/questions/1379564/c-tiny-memory-leak-with-stdmap
}

void GFBookkeeping::bookGFDetPlanes(const GFBKKey& key){
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<GFBKKey, std::vector<GFDetPlane> >::const_iterator it;
  it = fPlanes.find(key);
  if(it != fPlanes.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookGFDetPlanes()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fPlanes[key].resize(fNhits);// valgrind detects a memory leak here, but it is probably only due to stl optimizations: http://stackoverflow.com/questions/1379564/c-tiny-memory-leak-with-stdmap
}

void GFBookkeeping::bookNumbers(const GFBKKey& key, double val){ //val is set to 0. per default
  if(fNhits<0){
    GFException exc("fNhits not defined",__LINE__,__FILE__);
    throw exc;
  }
  std::map<GFBKKey, std::vector<double> >::const_iterator it;
  it = fNumbers.find(key);
  if(it != fNumbers.end()){
    std::ostringstream ostr;
    ostr << "The key " << key 
	 << " is already occupied in GFBookkeeping::bookNumbers()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  fNumbers[key].assign(fNhits, val);// valgrind detects a memory leak here, but it is probably only due to stl optimizations: http://stackoverflow.com/questions/1379564/c-tiny-memory-leak-with-stdmap
}


void GFBookkeeping::setVector(const GFBKKey& key, unsigned int index,
			      const TVectorD& vec){
  std::map<GFBKKey, std::vector<TVectorD> >::const_iterator it;
  it = fVectors.find(key);
  if(it == fVectors.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setVector()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setVector()";
   GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fVectors[key])[index].ResizeTo(vec);
  (fVectors[key])[index] = vec;
}

void GFBookkeeping::setMatrix(const GFBKKey& key, unsigned int index,
			    const TMatrixD& mat){
  std::map<GFBKKey, std::vector<TMatrixD> >::const_iterator it;
  it = fMatrices.find(key);
  if(it == fMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setMatrix()";
   GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fMatrices[key])[index].ResizeTo(mat);
  (fMatrices[key])[index] = mat;
}

void GFBookkeeping::setSymMatrix(const GFBKKey& key, unsigned int index,
				 const TMatrixDSym& mat){
  std::map<GFBKKey, std::vector<TMatrixDSym> >::const_iterator it;
  it = fSymMatrices.find(key);
  if(it == fSymMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setSymMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setSymMatrix()";
   GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fSymMatrices[key])[index].ResizeTo(mat);
  (fSymMatrices[key])[index] = mat;
}

void GFBookkeeping::setDetPlane(const GFBKKey& key, unsigned int index,
				  const GFDetPlane& pl){
  std::map<GFBKKey, std::vector<GFDetPlane> >::const_iterator it;
  it = fPlanes.find(key);
  if(it == fPlanes.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fPlanes[key])[index] = pl;
}

void GFBookkeeping::setNumber(const GFBKKey& key, unsigned int index,
			      const double& num){
  if(fNumbers[key].size() == 0){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::setNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::setNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  (fNumbers[key])[index] = num;
}


const TVectorD&
GFBookkeeping::getVector(const GFBKKey& key, unsigned int index) const {
  std::map<GFBKKey, std::vector<TVectorD> >::const_iterator it;
  it = fVectors.find(key);
  if(it == fVectors.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getVector()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getVector()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  return ((*it).second)[index];
}

const TMatrixD&
GFBookkeeping::getMatrix(const GFBKKey& key, unsigned int index) const {
  std::map<GFBKKey, std::vector<TMatrixD> >::const_iterator it;
  it = fMatrices.find(key);
  if(it == fMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  return ((*it).second)[index];
}

const TMatrixDSym&
GFBookkeeping::getSymMatrix(const GFBKKey& key, unsigned int index) const {
  std::map<GFBKKey, std::vector<TMatrixDSym> >::const_iterator it;
  it = fSymMatrices.find(key);
  if(it == fSymMatrices.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getSymMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getSymMatrix()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  return ((*it).second)[index];
}

const GFDetPlane&
GFBookkeeping::getDetPlane(const GFBKKey& key, unsigned int index) const {
  std::map<GFBKKey, std::vector<GFDetPlane> >::const_iterator it;
  it = fPlanes.find(key);
  if(it == fPlanes.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getGFDetPlane()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  return ((*it).second)[index];
}

double
GFBookkeeping::getNumber(const GFBKKey& key, unsigned int index) const {
  std::map<GFBKKey, std::vector<double> >::const_iterator it;
  it = fNumbers.find(key);
  if(it == fNumbers.end()){
    std::ostringstream ostr;
    ostr << "The key " << key << " is unknown in GFBookkeeping::getNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;
  }
  if(index>=(unsigned int)fNhits){
    std::ostringstream ostr;
    ostr << "The index " << index
	 << " is out of range in GFBookkeeping::getNumber()";
    GFException exc(ostr.str(),__LINE__,__FILE__);
    throw exc;    
  }
  return ((*it).second)[index];
}


std::vector< GFBKKey >
GFBookkeeping::getVectorKeys() const {
  std::vector< GFBKKey > keys;
  keys.reserve(fVectors.size());
  std::map<GFBKKey, std::vector<TVectorD> >::const_iterator it;
  for(it=fVectors.begin();it!=fVectors.end();++it){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}

std::vector< GFBKKey >
GFBookkeeping::getMatrixKeys() const {
  std::vector< GFBKKey > keys;
  keys.reserve(fMatrices.size());
  std::map<GFBKKey, std::vector<TMatrixD> >::const_iterator it;
  for(it=fMatrices.begin();it!=fMatrices.end();++it){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}

std::vector< GFBKKey >
GFBookkeeping::getSymMatrixKeys() const {
  std::vector< GFBKKey > keys;
  keys.reserve(fSymMatrices.size());
  std::map<GFBKKey, std::vector<TMatrixDSym> >::const_iterator it;
  for(it=fSymMatrices.begin();it!=fSymMatrices.end();++it){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}

std::vector< GFBKKey >
GFBookkeeping::getGFDetPlaneKeys() const {
  std::vector< GFBKKey > keys;
  keys.reserve(fPlanes.size());
  std::map<GFBKKey, std::vector<GFDetPlane> >::const_iterator it;
  for(it=fPlanes.begin();it!=fPlanes.end();++it){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}

std::vector< GFBKKey >
GFBookkeeping::getNumberKeys() const {
  std::vector< GFBKKey > keys;
  keys.reserve(fNumbers.size());
  std::map<GFBKKey, std::vector<double> >::const_iterator it;
  for(it=fNumbers.begin();it!=fNumbers.end();++it){
    if(it->second.size() != 0) keys.push_back(it->first);
  }
  return keys;
}


void GFBookkeeping::addFailedHit(unsigned int id){
  fFailedHits.push_back( id );
}

unsigned int GFBookkeeping::getNumFailed(){
  return fFailedHits.size();
}

unsigned int GFBookkeeping::hitFailed(unsigned int id){
  unsigned int retVal = 0;
  for(unsigned int i=0; i<fFailedHits.size(); ++i){
    if(fFailedHits[i] == id){
      ++retVal;
    }
  }
  return retVal;
}

void GFBookkeeping::clearFailedHits(){
  fFailedHits.clear();
}


void GFBookkeeping::reset() {

  clearFailedHits();

  unsigned int nHits(fNhits);
  if (fNhits < 0) nHits = 0; // needed because of std::vector::assign

  std::map<GFBKKey, std::vector<TVectorD > >::iterator itVec;
  for(itVec=fVectors.begin(); itVec!=fVectors.end(); ++itVec){
    if (itVec->second.size()==nHits) {
      for (unsigned int i=0; i<itVec->second.size(); ++i){
        itVec->second[i].Zero();
      }
    }
    else if (itVec->second.size()>0) {
      TVectorD vec(itVec->second[0].GetNrows()); // TVector filled with 0
      itVec->second.assign(nHits, vec);
    }
    else itVec->second.resize(nHits);
  }

  std::map<GFBKKey, std::vector<TMatrixD> >::iterator itMat;
  for(itMat=fMatrices.begin(); itMat!=fMatrices.end(); ++itMat){
    if (itMat->second.size()==nHits) {
      for (unsigned int i=0; i<itMat->second.size(); ++i){
        itMat->second[i].Zero();
      }
    }
    else if (itMat->second.size()>0) {
      TMatrixD mat(itMat->second[0].GetNrows(), itMat->second[0].GetNcols()); // TMatrix filled with 0
      itMat->second.assign(nHits, mat);
    }
    else itMat->second.resize(nHits);
  }

  std::map<GFBKKey, std::vector<TMatrixDSym> >::iterator itMatSym;
  for(itMatSym=fSymMatrices.begin(); itMatSym!=fSymMatrices.end(); ++itMatSym){
    if (itMatSym->second.size()==nHits) {
      for (unsigned int i=0; i<itMatSym->second.size(); ++i){
        itMatSym->second[i].Zero();
      }
    }
    else if (itMatSym->second.size()>0){
      TMatrixDSym matSym(itMatSym->second[0].GetNrows()); // TMatrix filled with 0
      itMatSym->second.assign(nHits, matSym);
    }
    else itMatSym->second.resize(nHits);
  }

  std::map<GFBKKey, std::vector<GFDetPlane> >::iterator itPl;
  for(itPl=fPlanes.begin(); itPl!=fPlanes.end(); ++itPl){
    if (itPl->second.size()==nHits) {
      for (unsigned int i=0; i<itPl->second.size(); ++i){
        itPl->second[i].reset();
      }
    }
    else {
      GFDetPlane pl;
      itPl->second.assign(nHits, pl);
    }
  }

  std::map<GFBKKey, std::vector<double> >::iterator itNum;
  for(itNum=fNumbers.begin(); itNum!=fNumbers.end(); ++itNum){
    itNum->second.assign(nHits, 0);
  }

}

void GFBookkeeping::clearAll(){
  fVectors.clear();
  fMatrices.clear();
  fSymMatrices.clear();
  fPlanes.clear();
  fNumbers.clear();
}


void GFBookkeeping::Print(const Option_t* option) const {
  std::cout << "=============GFBookkeeping::print()==============\n";
  std::cout << "GFBookkeeping has " << fNhits << " hits.\n";
  std::cout << "-----printing all vectors:------\n";
  std::vector<GFBKKey> keys = getVectorKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys[i] << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      getVector(keys[i],j).Print(option);
    }
  }
  std::cout << "-----printing all matrices:------\n";
  keys = getMatrixKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys[i] << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      getMatrix(keys[i],j).Print(option);
    }
  }
  std::cout << "-----printing all symmetric matrices:------\n";
  keys = getSymMatrixKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys[i] << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      getSymMatrix(keys[i],j).Print(option);
    }
  }
  std::cout << "-----printing all GFDetPlanes:------\n";
  keys = getGFDetPlaneKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys[i] << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      getDetPlane(keys[i],j).Print(option);
    }
  }
  std::cout << "-----printing all numbers:------\n";
  keys = getNumberKeys();
  for(unsigned int i=0;i<keys.size();++i){
    std::cout << "key " << keys[i] << " has " << fNhits
	      << " entries:\n";
    for(int j=0;j<fNhits;++j){
      std::cout << getNumber(keys[i],j) << std::endl;
    }
  }
  std::cout << "-----failed hits:------\n";
  for(unsigned int i=0;i<fFailedHits.size();++i){
    std::cout << fFailedHits[i] << " ";
  }
  std::cout << "==========GFBookkeeping::print() - Done==========" << std::endl;
}


ClassImp(GFBookkeeping)
