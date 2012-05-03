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
	fBlow = 1;
	fHitCovUpd = false;
	fHitCoordUpd = false;

}

GFAbsRecoHit* GFDafHit::getHit(unsigned int ihit) {

	return fRawHits.at(ihit);

}

void GFDafHit::setBlowUp(double blow_up) {

	fBlow = blow_up;
	fHitCovUpd = false;
	fHitCoordUpd = false;

}

void GFDafHit::setWeights(std::vector<double> weights) {

	fWeights = weights;
	fHitCovUpd = false;
	fHitCoordUpd = false;

}

const GFDetPlane& GFDafHit::getDetPlane(GFAbsTrackRep* rep) {

	return fRawHits.at(0)->getDetPlane(rep);

}

TMatrixT<double> GFDafHit::getHitCoord(const GFDetPlane& pl) {

	if(fHitCoordUpd && fPl == pl) return fHitCoord;

	if(fRawHits.size() == 1) {

		if(fHitCoord.GetNrows() ==0) fHitCoord.ResizeTo(fRawHits.at(0)->getHitCoord(pl));
		fHitCoord = fRawHits.at(0)->getHitCoord(pl);

	} else {

		if(fHitCovUpd != true || fPl != pl) getHitCov(pl);
		fHitCoord.Zero();
		for(unsigned int i=0;i<fRawHits.size();i++) {

			fHitCoord += fWeights.at(i) * fCovInvs.at(i) * fRawHits.at(i)->getHitCoord(pl);

		}

		fHitCoord = fHitCov * fHitCoord;
	}

	fPl = pl;
	fHitCoordUpd = true;
	return fHitCoord;

}

TMatrixT<double> GFDafHit::getHitCov(const GFDetPlane& pl) {

	if(fHitCovUpd && fPl == pl) return fHitCov;

	if(fHitCov.GetNrows() == 0) {
		fHitCov.ResizeTo(fRawHits.at(0)->getHitCov(pl));
		fHitCoord.ResizeTo(fRawHits.at(0)->getHitCoord(pl));
	}

	if(fRawHits.size() == 1) {

		if(((1/fWeights.at(0)) * fBlow) < pow(10,10)) {
			fHitCov = (1 / fWeights.at(0)) * fBlow * fRawHits.at(0)->getHitCov(pl);
		} else {
			fHitCov = pow(10,10) * fRawHits.at(0)->getHitCov(pl);
		}

	} else {

		fHitCov.Zero();
		fCovInvs.clear();

		for(unsigned int i=0;i<fRawHits.size();i++) {

			TMatrixT<double> CovInv;
			GFTools::invertMatrix((fBlow * fRawHits.at(i)->getHitCov(pl)), CovInv);
			fCovInvs.push_back(CovInv);
			fHitCov += fWeights.at(i) * CovInv;
		}

		TMatrixT<double> HitCovTemp(fHitCov);
		GFTools::invertMatrix(HitCovTemp, fHitCov);

	}

	fPl = pl;
	fHitCovUpd = true;
	return fHitCov;

}

TMatrixT<double> GFDafHit::getHMatrix(const GFAbsTrackRep* rep) {

	return fRawHits.at(0)->getHMatrix(rep);

}

GFDafHit* GFDafHit::clone() {

	GFDafHit* retval = new GFDafHit(fRawHits);
	retval->setWeights(fWeights);
	retval->setBlowUp(fBlow);
	return retval;

}

const std::string& GFDafHit::getPolicyName(){

	return fRawHits.at(0)->getPolicyName();

}

ClassImp(GFDafHit)
