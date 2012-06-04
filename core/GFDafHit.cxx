/* Copyright 2011, Technische Universitaet Muenchen,
Authors: Karl Bicker

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

#include<GFDafHit.h>

GFDafHit::GFDafHit(std::vector<GFAbsRecoHit*> HitsInPlane) {

	fRawHits = HitsInPlane;
	fWeights.assign(fRawHits.size(),1.);
	fHitUpd = false;

}

GFAbsRecoHit* GFDafHit::getHit(unsigned int ihit) {

	return fRawHits.at(ihit);

}

void GFDafHit::setWeights(std::vector<double> weights) {

	fWeights = weights;
	fHitUpd = false;
}

const GFDetPlane& GFDafHit::getDetPlane(GFAbsTrackRep* rep) {

	return fRawHits.at(0)->getDetPlane(rep);

}

void GFDafHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TMatrixT<double>& statePred,const TMatrixT<double>& covPred,TMatrixT<double>& m, TMatrixT<double>& V) {

  /*
  if(fHitUpd && fPl != pl) {
    GFException exc("GFDafHit::getMeasurement(): pl!=fPl",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  */

  if(fHitUpd && fPl == pl) {
    m.ResizeTo(fHitCoord);
    V.ResizeTo(fHitCov);
    m = fHitCoord;
    V = fHitCov;
    return;
  }

  if(fRawHits.size() == 1) {
    fRawHits.at(0)->getMeasurement(rep,pl,statePred,covPred,fHitCoord,fHitCov);
    static const double maxCovSize = 1.e10;
    if( 1.0/fWeights.at(0)  < maxCovSize) {
      fHitCov = (1.0 / fWeights.at(0)) * fHitCov;
    }
    else {
      fHitCov = maxCovSize * fHitCov;
    }
  } 

  else {
    //set the weighted-mean cov
    // this might seem like kind of a waste, but we need to make sure that fHitCov has the right dimensionality
    // and we dont know it from elsewhere
    fRawHits.at(0)->getMeasurement(rep,pl,statePred,covPred,fHitCoord,fHitCov);
    fHitCoord.Zero();
    fHitCov.Zero();
    fCovInvs.clear();
    TMatrixT<double> CovInv;
    TMatrixT<double>* coordTemp;
    TMatrixT<double> covTemp;
    std::vector<TMatrixT<double>* > coords;
    for(unsigned int i=0;i<fRawHits.size();i++) {
      coordTemp = new TMatrixT<double>;
      try{
	fRawHits.at(i)->getMeasurement(rep,pl,statePred,covPred,*coordTemp,covTemp);
	coords.push_back(coordTemp);
	GFTools::invertMatrix(covTemp, CovInv);
      }
      catch(GFException& e){
	for(unsigned int j=0;j<coords.size();++j) delete coords[j];
	delete coordTemp;
	throw e;
      }
      fCovInvs.push_back(CovInv);
      fHitCov += fWeights.at(i) * CovInv;
    }
    TMatrixT<double> HitCovTemp(fHitCov);

    try{
      GFTools::invertMatrix(HitCovTemp, fHitCov);
    }
    catch(GFException& e){
      for(unsigned int j=0;j<coords.size();++j) delete coords[j];
      throw e;
    }
    

    //set the weighted-mean coord
    //fRawHits.size()==coords.size() is a certainty here
    for(unsigned int i=0;i<fRawHits.size();i++) {
      fHitCoord += fWeights.at(i) * fCovInvs.at(i) * (*(coords.at(i)));
    }
    for(unsigned int j=0;j<coords.size();++j) delete coords[j];
    fHitCoord = fHitCov * fHitCoord;
  }
  //return by refernce
  m.ResizeTo(fHitCoord);
  V.ResizeTo(fHitCov);
  m = fHitCoord;
  V = fHitCov;
  fPl = pl;
  fHitUpd = true;
}


TMatrixT<double> GFDafHit::getHMatrix(const GFAbsTrackRep* rep) {

	return fRawHits.at(0)->getHMatrix(rep);

}

GFDafHit* GFDafHit::clone() {

	GFDafHit* retval = new GFDafHit(fRawHits);
	retval->setWeights(fWeights);
	return retval;

}

const std::string& GFDafHit::getPolicyName(){

	return fRawHits.at(0)->getPolicyName();

}

ClassImp(GFDafHit)
