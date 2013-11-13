/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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


#include "GFRaveVertex.h"
#include "GFRaveConverters.h"
#include <Exception.h>

#include <iostream>

namespace genfit {

//#define COUNT

#ifdef COUNT
static int instCount(0);
#endif

GFRaveVertex::GFRaveVertex() :
  cov_(3,3),
  ndf_(0),
  chi2_(0),
  id_(-1)
{
#ifdef COUNT
  std::cerr << "GFRaveVertex::GFRaveVertex() - Number of objects: " << ++instCount << std::endl;
#endif
}


GFRaveVertex::GFRaveVertex(const TVector3 & pos, const TMatrixDSym& cov,
                           const std::vector < GFRaveTrackParameters* > & smoothedTracks,
                           double ndf, double chi2, int id) :
  pos_(pos),
  cov_(cov),
  ndf_(ndf),
  chi2_(chi2),
  id_(id),
  smoothedTracks_(smoothedTracks)
{
  if (cov_.GetNrows()!=3 || cov_.GetNcols()!=3) {
    Exception exc("GFRaveVertex ==> Covariance is not 3x3!",__LINE__,__FILE__);
    throw exc;
  }

#ifdef COUNT
  std::cerr << "GFRaveVertex::GFRaveVertex(...) - Number of objects: " << ++instCount << std::endl;
#endif
}


GFRaveVertex::GFRaveVertex(const GFRaveVertex & vertex) :
  pos_(vertex.pos_),
  cov_(vertex.cov_),
  ndf_(vertex.ndf_),
  chi2_(vertex.chi2_),
  id_(vertex.id_)
{
  unsigned int nPar =  vertex.smoothedTracks_.size();
  smoothedTracks_.reserve(nPar);
  for (unsigned int i=0; i<nPar; ++i) {
    smoothedTracks_.push_back(new GFRaveTrackParameters(*(vertex.smoothedTracks_[i])));
  }

#ifdef COUNT
  std::cerr << "GFRaveVertex::GFRaveVertex(GFRaveVertex) - Number of objects: " << ++instCount << std::endl;
#endif
}


GFRaveVertex& GFRaveVertex::operator=(GFRaveVertex other) {
  swap(other);
  return *this;
}


void GFRaveVertex::swap(GFRaveVertex& other) {
  std::swap(this->pos_, other.pos_);
  this->cov_.ResizeTo(other.cov_);
  std::swap(this->cov_, other.cov_);
  std::swap(this->ndf_, other.ndf_);
  std::swap(this->chi2_, other.chi2_);
  std::swap(this->id_, other.id_);
  std::swap(this->smoothedTracks_, other.smoothedTracks_);
}


GFRaveVertex::~GFRaveVertex(){
  unsigned int nPar =  smoothedTracks_.size();
  for (unsigned int i=0; i<nPar; ++i) {
    delete smoothedTracks_[i];
  }

#ifdef COUNT
  std::cerr << "GFRaveVertex::~GFRaveVertex() - Number of objects: " << --instCount << std::endl;
#endif
}


void
GFRaveVertex::Print(const Option_t*) const {
  std::cout << "GFRaveVertex\n";
  std::cout << "Position: "; getPos().Print();
  std::cout << "Covariance: "; getCov().Print();
  std::cout << "Ndf: " << getNdf() << ", Chi2: " << getChi2() << ", Id: " << getId() << "\n";
  std::cout << "Number of tracks: " << getNTracks() << "\n";
  for (unsigned int i=0;  i<getNTracks(); ++i) {
    std::cout << " track " << i << ":\n"; getParameters(i)->Print();
  }
}

} /* End of namespace genfit */

