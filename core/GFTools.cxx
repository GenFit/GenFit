
#include "GFTools.h"

#include <cmath>
#include <memory>
#include <typeinfo>

#include <TDecompChol.h>
#include <TMath.h>

TVectorD GFTools::getSmoothedPos(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TVectorD smoothed_state;
	TMatrixDSym smoothed_cov;
	TVectorD pos;

	if(GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov)) {

		const TMatrixD& H( trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep)) );
		TVectorD pos_tmp(H * smoothed_state);
		pos.ResizeTo(pos_tmp);
		pos = pos_tmp;

	}
	return pos;

}


TVector3 GFTools::getSmoothedPosXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorD smoothed_state;
  TMatrixDSym smoothed_cov;
  TMatrixD auxInfo;
  GFDetPlane smoothing_plane;

  getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getPos(smoothing_plane);
}

TVector3 GFTools::getBiasedSmoothedPosXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorD smoothed_state;
  TMatrixDSym smoothed_cov;
  TMatrixD auxInfo;
  GFDetPlane smoothing_plane;

  getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getPos(smoothing_plane);
}


TVector3 GFTools::getSmoothedMomXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorD smoothed_state;
  TMatrixDSym smoothed_cov;
  TMatrixD auxInfo;
  GFDetPlane smoothing_plane;

  getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getMom(smoothing_plane);
}

TVector3 GFTools::getBiasedSmoothedMomXYZ(const GFTrack* trk, unsigned int irep, unsigned int ihit){

  std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

  TVectorD smoothed_state;
  TMatrixDSym smoothed_cov;
  TMatrixD auxInfo;
  GFDetPlane smoothing_plane;

  getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

  if(rep->hasAuxInfo()) {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
  } else {
    rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
  }

  return rep->getMom(smoothing_plane);
}


TVectorD GFTools::getBiasedSmoothedPos(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

  TVectorD smoothed_state;
  TMatrixDSym smoothed_cov;
  TVectorD pos;

	if(GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov)) {

	  const TMatrixD& H( trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep)) );
		//H.Print();smoothed_state.Print();
		TVectorD pos_tmp(H * smoothed_state);
		pos.ResizeTo(pos_tmp);
		pos = pos_tmp;

	}
	return pos;

}

TMatrixDSym GFTools::getSmoothedCov(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TVectorD smoothed_state;
	TMatrixDSym smoothed_cov;
	GFDetPlane pl;

	GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl);
	const TMatrixD& H( trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep)) );

	TMatrixDSym cov(smoothed_cov);
	cov.Similarity(H);

	return cov;

}

TMatrixDSym GFTools::getBiasedSmoothedCov(const GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TVectorD smoothed_state;
	TMatrixDSym smoothed_cov;
	GFDetPlane pl;

	GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl);
	const TMatrixD& H( trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep)) );
	TMatrixDSym cov(smoothed_cov);
	cov.Similarity(H);

	return cov;

}

bool GFTools::getSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorD& smoothed_state, TMatrixDSym& smoothed_cov) {

	GFDetPlane pl;
	TMatrixD auxInfo;
	return GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl, auxInfo);

}

bool GFTools::getBiasedSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorD& smoothed_state, TMatrixDSym& smoothed_cov) {

	GFDetPlane pl;
	TMatrixD auxInfo;
	return GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl, auxInfo);

}

bool GFTools::getSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorD& smoothed_state, TMatrixDSym& smoothed_cov, GFDetPlane& smoothing_plane) {

	TMatrixD auxInfo;
	return GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

}

bool GFTools::getBiasedSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorD& smoothed_state, TMatrixDSym& smoothed_cov, GFDetPlane& smoothing_plane) {

	TMatrixD auxInfo;
	return GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

}

bool GFTools::getSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorD& smoothed_state, TMatrixDSym& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixD& auxInfo) {

	if(!trk->getSmoothing()) {
		std::cout << "Trying to get smoothed hit coordinates from a track without smoothing! Aborting..." << std::endl;
		return false;
	}

	if(ihit >= trk->getNumHits()) {
		std::cout << "Hit number out of bounds while trying to get smoothed hit coordinates! Aborting..." <<std::endl;
		return false;
	}

	unsigned int dim(trk->getTrackRep(irep)->getDim());
	smoothed_state.ResizeTo(dim);
	smoothed_cov.ResizeTo(dim, dim);

	// get aux info
	if(trk->getTrackRep(irep)->hasAuxInfo()) {
    auxInfo.ResizeTo(trk->getBK(irep)->getMatrix(GFBKKey_fAuxInfo,ihit));
    auxInfo = trk->getBK(irep)->getMatrix(GFBKKey_fAuxInfo,ihit);
	}

	if(ihit == 0) { // first hit -> get prediction from backward filter
	  smoothed_state = trk->getBK(irep)->getVector(GFBKKey_bSt,ihit);
	  smoothed_cov = trk->getBK(irep)->getSymMatrix(GFBKKey_bCov,ihit);
	  smoothing_plane = trk->getBK(irep)->getDetPlane(GFBKKey_bPl,ihit);
		return true;
	}

	if(ihit == trk->getNumHits()-1) { // last hit -> get prediction from forward filter
	  smoothed_state = trk->getBK(irep)->getVector(GFBKKey_fSt,ihit);
	  smoothed_cov = trk->getBK(irep)->getSymMatrix(GFBKKey_fCov,ihit);
	  smoothing_plane = trk->getBK(irep)->getDetPlane(GFBKKey_fPl,ihit);
		return true;
	}


	// calculate "mean" between forward and backward predictions
  const TVectorD& fSt(trk->getBK(irep)->getVector(GFBKKey_fSt, ihit));
  const TMatrixDSym& fCov(trk->getBK(irep)->getSymMatrix(GFBKKey_fCov, ihit));
  smoothing_plane = trk->getBK(irep)->getDetPlane(GFBKKey_fPl, ihit);

  TVectorD bSt(trk->getBK(irep)->getVector(GFBKKey_bSt, ihit)); // can't be const ref since it possibly has to be altered in extrapolate
  TMatrixDSym bCov(trk->getBK(irep)->getSymMatrix(GFBKKey_bCov, ihit)); // can't be const ref since it possibly has to be altered in extrapolate
  const GFDetPlane& bPl(trk->getBK(irep)->getDetPlane(GFBKKey_bPl, ihit));

  if(smoothing_plane != bPl) {
    // if the two planes are not identical, we actually would need to extrapolate back
    // from the update of ihit+1.
    // But since the planes should be close, I think just extrapolating
    // from the prediction at ihit to the plane is a good approximation.
    // Eventually, if a backwards propagation is done, the covariance can be a bit too large.
    // But it should be much faster this way.

    const TMatrixD* bAuxInfoP;

    std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

    if(trk->getTrackRep(irep)->hasAuxInfo()) {
      bAuxInfoP = &(trk->getBK(irep)->getMatrix(GFBKKey_bAuxInfo, ihit));
    } else {
      bAuxInfoP = NULL;
    }

    rep->setData(bSt, bPl, &bCov, bAuxInfoP);
    rep->extrapolate(smoothing_plane, bSt, bCov);
  }

	TMatrixDSym fCovInvert;
	TMatrixDSym bCovInvert;

	GFTools::invertMatrix(fCov, fCovInvert);
	GFTools::invertMatrix(bCov, bCovInvert);

	GFTools::invertMatrix(fCovInvert + bCovInvert, smoothed_cov);

	smoothed_state = smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt);

	return true;

}

bool GFTools::getBiasedSmoothedData(const GFTrack* trk, unsigned int irep, unsigned int ihit, TVectorD& smoothed_state, TMatrixDSym& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixD& auxInfo) {

	if(!trk->getSmoothing()) {
		std::cout << "Trying to get smoothed hit coordinates from a track without smoothing! Aborting..." << std::endl;
		TMatrixD failed(1,1);
		return false;
	}

	if(ihit >= trk->getNumHits()) {
		std::cout << "Hit number out of bounds while trying to get smoothed hit coordinates! Aborting..." <<std::endl;
		TMatrixD failed(1,1);
		failed(0,0) = 0;
		return false;
	}

	unsigned int dim(trk->getTrackRep(irep)->getDim());
	smoothed_state.ResizeTo(dim);
	smoothed_cov.ResizeTo(dim, dim);


	if(trk->getTrackRep(irep)->hasAuxInfo()) {
	  auxInfo.ResizeTo(trk->getBK(irep)->getMatrix(GFBKKey_fAuxInfo,ihit));
	  auxInfo = trk->getBK(irep)->getMatrix(GFBKKey_fAuxInfo,ihit);
	}

	if(ihit == 0) { // first hit -> get update from backward filter
	  smoothed_state = trk->getBK(irep)->getVector(GFBKKey_bUpSt,ihit);
	  smoothed_cov = trk->getBK(irep)->getSymMatrix(GFBKKey_bUpCov,ihit);
	  smoothing_plane = trk->getBK(irep)->getDetPlane(GFBKKey_bPl,ihit);
		return true;
	}

	if(ihit == trk->getNumHits()-1) { // last hit -> get update from forward filter
	  smoothed_state = trk->getBK(irep)->getVector(GFBKKey_fUpSt,ihit);
	  smoothed_cov = trk->getBK(irep)->getSymMatrix(GFBKKey_fUpCov,ihit);
	  smoothing_plane = trk->getBK(irep)->getDetPlane(GFBKKey_fPl,ihit);
		return true;
	}


	// calculate "mean" between forward update and backward prediction
	const TVectorD& fUpSt(trk->getBK(irep)->getVector(GFBKKey_fUpSt,ihit));
	const TMatrixDSym& fUpCov(trk->getBK(irep)->getSymMatrix(GFBKKey_fUpCov,ihit));
	smoothing_plane = trk->getBK(irep)->getDetPlane(GFBKKey_fPl,ihit);

	TVectorD bSt(trk->getBK(irep)->getVector(GFBKKey_bSt,ihit)); // can't be const ref since it possibly has to be altered in extrapolate
	TMatrixDSym bCov(trk->getBK(irep)->getSymMatrix(GFBKKey_bCov,ihit)); // can't be const ref since it possibly has to be altered in extrapolate
	const GFDetPlane& bPl(trk->getBK(irep)->getDetPlane(GFBKKey_bPl,ihit));


	if(smoothing_plane != bPl) {
    // if the two planes are not identical, we actually would need to extrapolate back
    // from the update of ihit+1.
    // But since the planes should be close, I think just extrapolating
    // from the prediction at ihit to the plane is a good approximation.
    // Eventually, if a backwards propagation is done, the covariance can be a bit too large.
    // But it should be much faster this way.

    const TMatrixD* bAuxInfoP;

    std::auto_ptr<GFAbsTrackRep> rep(trk->getTrackRep(irep)->clone());

    if(trk->getTrackRep(irep)->hasAuxInfo()) {
      bAuxInfoP = &(trk->getBK(irep)->getMatrix(GFBKKey_bAuxInfo, ihit));
    } else {
      bAuxInfoP = NULL;
    }

    rep->setData(bSt, bPl, &bCov, bAuxInfoP);
    rep->extrapolate(smoothing_plane, bSt, bCov);

	}

	TMatrixDSym fCovInvert;
	TMatrixDSym bCovInvert;

	GFTools::invertMatrix(fUpCov, fCovInvert);
	GFTools::invertMatrix(bCov, bCovInvert);

	GFTools::invertMatrix(fCovInvert + bCovInvert, smoothed_cov);

	smoothed_state = smoothed_cov * (fCovInvert*fUpSt + bCovInvert*bSt);

	return true;

}

const GFDetPlane& GFTools::getSmoothingPlane(const GFTrack* trk, unsigned int irep, unsigned int ihit) {
	return trk->getBK(irep)->getDetPlane(GFBKKey_fPl,ihit);
}

double GFTools::getTrackLength(const GFTrack* trk, unsigned int irep, unsigned int startHit, unsigned int endHit){

  if(!trk->getSmoothing()) {
    GFException exc("Trying to get tracklength from a track without smoothing!",__LINE__,__FILE__);
    throw exc;
  }

  if(startHit >= trk->getNumHits() || endHit >= trk->getNumHits()) {
    GFException exc("Hit number out of bounds while trying to get tracklength!",__LINE__,__FILE__);
    throw exc;
  }

  bool inv(false);
  if(startHit > endHit) {
    unsigned int biggerOne = startHit;
    startHit = endHit;
    endHit = biggerOne;
    inv = true;
  }

  if (startHit==0 && endHit==0) endHit = trk->getNumHits()-1;
  if (startHit == endHit) return 0.;

  double totLen(0), fLen(0), bLen(0);

  for (unsigned int i=startHit; i!=endHit; ++i){
    fLen = trk->getBK(irep)->getNumber(GFBKKey_fExtLen, i+1);
    bLen = trk->getBK(irep)->getNumber(GFBKKey_bExtLen, i);
    totLen += 0.5*(fLen - bLen);
  }

  if (inv) return -totLen;
  return totLen;

}

void GFTools::invertMatrix(const TMatrixDSym& mat, TMatrixDSym& inv, double* determinant){
	inv.ResizeTo(mat);

	// check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
	if (!(mat<1.E100) || !(mat>-1.E100)){
		GFException e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1e100)",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;	
	}
	// do the trivial inversions for 1x1 and 2x2 matrices manually
	if (mat.GetNrows() == 1){
	  if (determinant != NULL) *determinant = mat(0,0);
	  inv(0,0) = 1./mat(0,0);
	  return;
	} if (mat.GetNrows() == 2){
	  double det = mat(0,0)*mat(1,1) - mat(1,0)*mat(1,0);
	  if (determinant != NULL) *determinant = det;
	  if(fabs(det) < 1E-50){
	    GFException e("cannot invert matrix GFTools::invertMatrix(), determinant = 0",
	        __LINE__,__FILE__);
	    e.setFatal();
	    throw e;
	  }
	  det = 1./det;
	  inv(0,0) =             det * mat(1,1);
	  inv(0,1) = inv(1,0) = -det * mat(1,0);
	  inv(1,1) =             det * mat(0,0);
	  return;
	}

	// else use TDecompChol
	bool status = 0;
	TDecompChol invertAlgo(mat, 1E-50);

	status = invertAlgo.Invert(inv);
	if(status == 0){
		GFException e("cannot invert matrix GFTools::invertMatrix(), status = 0",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}

	if (determinant != NULL) {
    double d1, d2;
    invertAlgo.Det(d1, d2);
	  *determinant = ldexp(d1, d2);
	}
}

void GFTools::invertMatrix(TMatrixDSym& mat, double* determinant){
	// check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
	if (!(mat<1.E100) || !(mat>-1.E100)){
		GFException e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1e100)",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;	
	}
	// do the trivial inversions for 1x1 and 2x2 matrices manually
	if (mat.GetNrows() == 1){
	  if (determinant != NULL) *determinant = mat(0,0);
	  mat(0,0) = 1./mat(0,0);
	  return;
	} else if (mat.GetNrows() == 2){
	  double *arr = mat.GetMatrixArray();
	  double det = arr[0]*arr[3] - arr[1]*arr[1];
	  if (determinant != NULL) *determinant = det;
	  if(fabs(det) < 1E-50){
	    GFException e("cannot invert matrix GFTools::invertMatrix(), determinant = 0",
	        __LINE__,__FILE__);
	    e.setFatal();
	    throw e;
	  }
	  det = 1./det;
	  double temp[3];
	  temp[0] =  det * arr[3];
	  temp[1] = -det * arr[1];
	  temp[2] =  det * arr[0];
	  //double *arr = mat.GetMatrixArray();
	  arr[0] = temp[0];
	  arr[1] = arr[2] = temp[1];
	  arr[3] = temp[2];
	  return;
	}

	// else use TDecompChol
	bool status = 0;
	TDecompChol invertAlgo(mat, 1E-50);

	status = invertAlgo.Invert(mat);
	if(status == 0){
		GFException e("cannot invert matrix GFTools::invertMatrix(), status = 0",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}

  if (determinant != NULL) {
    double d1, d2;
    invertAlgo.Det(d1, d2);
    *determinant = ldexp(d1, d2);
  }
}

void GFTools::updateRepSmoothed(GFTrack* trk, unsigned int irep, unsigned int ihit) {

  TVectorD smoothed_state;
  TMatrixD auxInfo;
	TMatrixDSym smoothed_cov;
	GFDetPlane smoothing_plane;

	getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

	GFAbsTrackRep* rep = trk->getTrackRep(irep);
	if(rep->hasAuxInfo()) {
		rep->setData(smoothed_state, smoothing_plane, &smoothed_cov, &auxInfo);
	} else {
		rep->setData(smoothed_state, smoothing_plane, &smoothed_cov);
	}

	trk->setRepAtHit(irep, ihit);
}

double GFTools::getSmoothedChiSqu(const GFTrack* trk, unsigned int irep, unsigned int ihit){
	TVectorD smoothed_state;
	TMatrixDSym smoothed_cov;
	GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov);
	const TMatrixD& H( trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep)) );
	TVectorD pos = H * smoothed_state;
	const TVectorD& m( trk->getHit(ihit)->getRawHitCoord() ); //measurement of hit
	const TMatrixDSym& V( trk->getHit(ihit)->getRawHitCov() ); //covariance matrix of hit
	TVectorD res = m - pos;
	TMatrixDSym R = V - smoothed_cov.Similarity(H);
	TMatrixDSym invR;
	invertMatrix(R,invR);
	double smoothedChiSqu = invR.Similarity(res);
	return smoothedChiSqu;
}

unsigned int GFTools::getClosestHit(const GFTrack* trk, unsigned int irep, const TVector3& pos, double& distance, bool checkEveryHit){

  if(!trk->getSmoothing()) {
    GFException exc("Trying to get closest hit from a track without smoothing!",__LINE__,__FILE__);
    throw exc;
  }

	unsigned int nHits = trk->getNumHits();
	if (nHits == 1) return 0;

	TVector3 hitPos;
	double minDist(1.E99), dist;
	unsigned int minId(0);

	if (checkEveryHit || nHits<4){
		// brute force search
		for (unsigned int i=0; i<nHits; ++i){
			hitPos = getSmoothedPosXYZ(trk, irep, i);
			dist = (pos-hitPos).Mag();
			if (dist<minDist) {
				minDist = dist;
				minId = i;
			}
		}
	}
	else { // hill climbing algorithm
		double distFirst = (pos-getSmoothedPosXYZ(trk, irep, 0)).Mag();
		double distLast = (pos-getSmoothedPosXYZ(trk, irep, nHits-1)).Mag();
		if (distFirst <= distLast){
			minDist = distFirst;
			minId = 0;
			for (unsigned int i=1; i<nHits-1; ++i){
				hitPos = getSmoothedPosXYZ(trk, irep, i);
				dist = (pos-hitPos).Mag();
				if (dist<=minDist) {
					minDist = dist;
					minId = i;
				}
				else break; // distance is increasing again!
			}
		}
		else {
			minDist = distLast;
			minId = nHits-1;
			for (unsigned int i=nHits-2; i>1; --i){
				hitPos = getSmoothedPosXYZ(trk, irep, i);
				dist = (pos-hitPos).Mag();
				if (dist<=minDist) {
					minDist = dist;
					minId = i;
				}
				else break; // distance is increasing again!
			}
		}
	}

	distance = minDist;
	return minId;
}

