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


#include "GFRaveVertex.h"
#include "GFRaveConverters.h"
#include <GFException.h>

#include <iostream>

using namespace std;


GFRaveVertex::GFRaveVertex() :
  fCov(3,3),
  fNdf(0),
  fChi2(0),
  fId(-1)
{
  fSmoothedTracks = new TObjArray();
}


GFRaveVertex::GFRaveVertex(TVector3 pos, TMatrixT<double> cov,
                           std::vector < GFRaveTrackParameters* > smoothedTracks,
                           double ndf, double chi2, int id) :
  fPos(pos),
  fCov(cov),
  fNdf(ndf),
  fChi2(chi2),
  fId(id)
{
  if (fCov.GetNrows()!=3 || fCov.GetNcols()!=3) {
    GFException exc("GFRaveVertex ==> Covariance is not 3x3!",__LINE__,__FILE__);
    throw exc;
  }

  unsigned int nTrks(smoothedTracks.size());
  fSmoothedTracks = new TObjArray(nTrks);
  for (unsigned int i=0; i<nTrks; ++i){
    fSmoothedTracks->Add(smoothedTracks[i]);
  }
}


GFRaveVertex::GFRaveVertex(const GFRaveVertex & vertex){
  fPos = vertex.fPos;
  fCov.ResizeTo(vertex.fCov);
  fCov = vertex.fCov;
  fNdf = vertex.fNdf;
  fChi2 = vertex.fChi2;
  fId = vertex.fId;
  fSmoothedTracks = (TObjArray*)(vertex.fSmoothedTracks->Clone());
}


GFRaveVertex& GFRaveVertex::operator=(const GFRaveVertex & vertex) {
  if(fSmoothedTracks!=NULL){
    fSmoothedTracks->Delete();
    delete fSmoothedTracks;
  }

  fPos = vertex.fPos;
  fCov.ResizeTo(vertex.fCov);
  fCov = vertex.fCov;
  fNdf = vertex.fNdf;
  fChi2 = vertex.fChi2;
  fId = vertex.fId;
  fSmoothedTracks = (TObjArray*)(vertex.fSmoothedTracks->Clone());
}


GFRaveVertex::~GFRaveVertex(){
  if(fSmoothedTracks!=NULL){
    fSmoothedTracks->Delete();
    delete fSmoothedTracks;
  }
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

