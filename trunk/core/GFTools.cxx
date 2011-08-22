
#include <GFTools.h>
#include <typeinfo>
TMatrixT<double> GFTools::getSmoothedPos(GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TMatrixT<double> smoothed_state;
	TMatrixT<double> smoothed_cov;
	TMatrixT<double> pos;

	if(GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov)) {

		TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
		TMatrixT<double> pos_tmp(H * smoothed_state);
		pos.ResizeTo(pos_tmp);
		pos = pos_tmp;

	}
	return pos;

}

TMatrixT<double> GFTools::getBiasedSmoothedPos(GFTrack* trk, unsigned int irep, unsigned int ihit) {

TMatrixT<double> smoothed_state;
	TMatrixT<double> smoothed_cov;
	TMatrixT<double> pos;

	if(GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov)) {

		TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
		H.Print();smoothed_state.Print();
		TMatrixT<double> pos_tmp(H * smoothed_state);
		pos.ResizeTo(pos_tmp);
		pos = pos_tmp;

	}
	return pos;

}

TMatrixT<double> GFTools::getSmoothedCov(GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TMatrixT<double> smoothed_state;
	TMatrixT<double> smoothed_cov;
	GFDetPlane pl;

	GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl);
	TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));

	TMatrixT<double> cov_tmp(smoothed_cov,TMatrixT<double>::kMultTranspose,H);
	TMatrixT<double> cov(H,TMatrixT<double>::kMult,cov_tmp);

	return cov;

}

TMatrixT<double> GFTools::getBiasedSmoothedCov(GFTrack* trk, unsigned int irep, unsigned int ihit) {

	TMatrixT<double> smoothed_state;
	TMatrixT<double> smoothed_cov;
	GFDetPlane pl;

	GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl);
	TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));

	TMatrixT<double> cov_tmp(smoothed_cov,TMatrixT<double>::kMultTranspose,H);
	TMatrixT<double> cov(H,TMatrixT<double>::kMult,cov_tmp);

	return cov;

}

bool GFTools::getSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov) {

	GFDetPlane pl;
	TMatrixT<double> auxInfo;
	return GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl, auxInfo);

}

bool GFTools::getBiasedSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov) {

	GFDetPlane pl;
	TMatrixT<double> auxInfo;
	return GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, pl, auxInfo);

}

bool GFTools::getSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane) {

	TMatrixT<double> auxInfo;
	return GFTools::getSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

}

bool GFTools::getBiasedSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane) {

	TMatrixT<double> auxInfo;
	return GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov, smoothing_plane, auxInfo);

}

bool GFTools::getSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixT<double>& auxInfo) {

	if(!trk->getSmoothing()) {
		std::cout << "Trying to get smoothed hit coordinates from a track without smoothing! Aborting..." << std::endl;
		TMatrixT<double> failed(1,1);
		return false;
	}

	if(ihit >= trk->getNumHits()) {
		std::cout << "Hit number out of bounds while trying to get smoothed hit coordinates! Aborting..." <<std::endl;
		TMatrixT<double> failed(1,1);
		failed(0,0) = 0;
		return false;
	}

	TMatrixT<double> fUpSt;
	TMatrixT<double> fUpCov;
	TMatrixT<double> fAuxInfo;
	TMatrixT<double>* fAuxInfoP;
	TMatrixT<double> bUpSt;
	TMatrixT<double> bUpCov;
	TMatrixT<double> bAuxInfo;
	TMatrixT<double>* bAuxInfoP;

	GFDetPlane fPl;
	GFDetPlane bPl;

	GFAbsTrackRep* rep = trk->getTrackRep(irep)->clone();

	if(trk->getTrackRep(irep)->hasAuxInfo()) {
		trk->getBK(irep)->getMatrix("fAuxInfo",ihit,auxInfo);
	}

	if(ihit == 0) {
		trk->getBK(irep)->getMatrix("bUpSt",ihit+1,bUpSt);
		trk->getBK(irep)->getMatrix("bUpCov",ihit+1,bUpCov);
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);
		if(trk->getTrackRep(irep)->hasAuxInfo()) {
			trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
			bAuxInfoP = &bAuxInfo;
		} else {
			bAuxInfoP = NULL;
		}
		if(bUpSt.GetNrows() == 0) return false;
		rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
		rep->extrapolate(smoothing_plane,smoothed_state,smoothed_cov);
		return true;
	}

	if(ihit == trk->getNumHits()-1) {
		trk->getBK(irep)->getMatrix("fUpSt",ihit-1,fUpSt);
		trk->getBK(irep)->getMatrix("fUpCov",ihit-1,fUpCov);
		trk->getBK(irep)->getDetPlane("fPl",ihit-1,fPl);
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		if(trk->getTrackRep(irep)->hasAuxInfo()) {
			trk->getBK(irep)->getMatrix("fAuxInfo",ihit-1,fAuxInfo);
			fAuxInfoP = &fAuxInfo;
		} else {
			fAuxInfoP = NULL;
		}
		if(fUpSt.GetNrows() == 0) return false;
		rep->setData(fUpSt,fPl,&fUpCov,fAuxInfoP);
		rep->extrapolate(smoothing_plane,smoothed_state,smoothed_cov);
		return true;
	}

	trk->getBK(irep)->getMatrix("fUpSt",ihit-1,fUpSt);
	trk->getBK(irep)->getMatrix("fUpCov",ihit-1,fUpCov);
	trk->getBK(irep)->getMatrix("bUpSt",ihit+1,bUpSt);
	trk->getBK(irep)->getMatrix("bUpCov",ihit+1,bUpCov);
	if(trk->getTrackRep(irep)->hasAuxInfo()) {
		trk->getBK(irep)->getMatrix("fAuxInfo",ihit-1,fAuxInfo);
		trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
		fAuxInfoP = &fAuxInfo;
		bAuxInfoP = &bAuxInfo;
	} else {
		fAuxInfoP = bAuxInfoP = NULL;
	}
	trk->getBK(irep)->getDetPlane("fPl",ihit-1,fPl);
	trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
	trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);


	TMatrixT<double> fSt;
	TMatrixT<double> fCov;
	TMatrixT<double> bSt;
	TMatrixT<double> bCov;

	if(fUpSt.GetNrows() == 0 || bUpSt.GetNrows() == 0) return false;

	rep->setData(fUpSt,fPl,&fUpCov,fAuxInfoP);
	rep->extrapolate(smoothing_plane,fSt,fCov);
	rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
	rep->extrapolate(smoothing_plane,bSt,bCov);

	delete rep;

	TMatrixT<double> fCovInvert;
	TMatrixT<double> bCovInvert;

	GFTools::invertMatrix(fCov, fCovInvert);
	GFTools::invertMatrix(bCov, bCovInvert);

	GFTools::invertMatrix(fCovInvert + bCovInvert, smoothed_cov);

	smoothed_state.ResizeTo(smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt));
	smoothed_state = smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt);

	return true;

}

bool GFTools::getBiasedSmoothedData(GFTrack* trk, unsigned int irep, unsigned int ihit, TMatrixT<double>& smoothed_state, TMatrixT<double>& smoothed_cov, GFDetPlane& smoothing_plane, TMatrixT<double>& auxInfo) {

	if(!trk->getSmoothing()) {
		std::cout << "Trying to get smoothed hit coordinates from a track without smoothing! Aborting..." << std::endl;
		TMatrixT<double> failed(1,1);
		return false;
	}

	if(ihit >= trk->getNumHits()) {
		std::cout << "Hit number out of bounds while trying to get smoothed hit coordinates! Aborting..." <<std::endl;
		TMatrixT<double> failed(1,1);
		failed(0,0) = 0;
		return false;
	}

	TMatrixT<double> bUpSt;
	TMatrixT<double> bUpCov;
	TMatrixT<double> bAuxInfo;
	TMatrixT<double>* bAuxInfoP;
	GFDetPlane bPl;

	GFAbsTrackRep* rep = trk->getTrackRep(irep)->clone();
	
	if(trk->getTrackRep(irep)->hasAuxInfo()) {
		trk->getBK(irep)->getMatrix("fAuxInfo",ihit,auxInfo);
	}

	if(ihit == 0) {
		trk->getBK(irep)->getMatrix("bUpSt",ihit,smoothed_state);
		trk->getBK(irep)->getMatrix("bUpCov",ihit,smoothed_cov);
		trk->getBK(irep)->getDetPlane("bPl",ihit,smoothing_plane);
		return true;
	}

	if(ihit == trk->getNumHits()-1) {
		trk->getBK(irep)->getMatrix("fUpSt",ihit,smoothed_state);
		trk->getBK(irep)->getMatrix("fUpCov",ihit,smoothed_cov);
		trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
		return true;
	}

	TMatrixT<double> fSt;
	TMatrixT<double> fCov;

	trk->getBK(irep)->getMatrix("fUpSt",ihit,fSt);
	trk->getBK(irep)->getMatrix("fUpCov",ihit,fCov);
	trk->getBK(irep)->getMatrix("bUpSt",ihit+1,bUpSt);
	trk->getBK(irep)->getMatrix("bUpCov",ihit+1,bUpCov);
	if(trk->getTrackRep(irep)->hasAuxInfo()) {
		trk->getBK(irep)->getMatrix("bAuxInfo",ihit+1,bAuxInfo);
		bAuxInfoP = &bAuxInfo;
	} else {
		bAuxInfoP = NULL;
	}
	trk->getBK(irep)->getDetPlane("fPl",ihit,smoothing_plane);
	trk->getBK(irep)->getDetPlane("bPl",ihit+1,bPl);


	TMatrixT<double> bSt;
	TMatrixT<double> bCov;

	if(bUpSt.GetNrows() == 0) return false;

	rep->setData(bUpSt,bPl,&bUpCov,bAuxInfoP);
	rep->extrapolate(smoothing_plane,bSt,bCov);

	delete rep;

	TMatrixT<double> fCovInvert;
	TMatrixT<double> bCovInvert;

	GFTools::invertMatrix(fCov, fCovInvert);
	GFTools::invertMatrix(bCov, bCovInvert);

	GFTools::invertMatrix(fCovInvert + bCovInvert, smoothed_cov);

	smoothed_state.ResizeTo(smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt));
	smoothed_state = smoothed_cov * (fCovInvert*fSt + bCovInvert*bSt);

	return true;

}

GFDetPlane GFTools::getSmoothingPlane(GFTrack* trk, unsigned int irep, unsigned int ihit) {

	GFDetPlane pl;
	trk->getBK(irep)->getDetPlane("fPl",ihit,pl);
	return pl;

}


/*
void GFTools::invertMatrix(const TMatrixT<double>& mat, TMatrixT<double>& inv){
	inv.ResizeTo(mat);
	// since ROOT has issues with absolute matrix entries that are below 10**-16 we force the 11 element to be one
	//double factor = 1./mat[0][0];
	//inv = (factor*mat);
	inv = mat;
	inv.SetTol(1E-150);
	double det=0.;
	inv.Invert(&det);
	if(TMath::IsNaN(det)) {
		GFException e("GFTools::invertMatrix(): det of matrix is nan",__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}
	if(det==0){
		GFException e("cannot invert matrix GFTools::invertMatrix() - det=0",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}
	//recover multiplication with factor at the beginning
	//inv *= factor;
}
*/

void GFTools::invertMatrix(const TMatrixT<double>& mat, TMatrixT<double>& inv){
	inv.ResizeTo(mat);
	bool status = 0;
	TDecompSVD invertAlgo(mat);

  // check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
	if (!(mat<1.E100) || !(mat>-1.E100)){
		GFException e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1E100)",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;	
	}
	
	invertAlgo.SetTol(1E-50); //this is a hack because a tolerance of 1E-22 does not make any sense for doubles only for floats
	
	inv = invertAlgo.Invert(status);
	if(status == 0){
		GFException e("cannot invert matrix GFTools::invertMatrix(), status = 0",
				__LINE__,__FILE__);
		e.setFatal();
		throw e;
	}
}

double GFTools::getSmoothedChiSqu(GFTrack* const trk, unsigned int irep, unsigned int ihit){
	TMatrixT<double> smoothed_state;
	TMatrixT<double> smoothed_cov;
	GFTools::getBiasedSmoothedData(trk, irep, ihit, smoothed_state, smoothed_cov);
	TMatrixT<double> H = trk->getHit(ihit)->getHMatrix(trk->getTrackRep(irep));
	TMatrixT<double> HT(TMatrixT<double>::kTransposed, H);
	TMatrixT<double> pos = H * smoothed_state;
	TMatrixT<double> m = trk->getHit(ihit)->getRawHitCoord(); //measurement of hit
	TMatrixT<double> V = trk->getHit(ihit)->getRawHitCov(); //covariance matrix of hit
	TMatrixT<double> res = m - pos;
	TMatrixT<double> resT(TMatrixT<double>::kTransposed, res);
	TMatrixT<double> R = V - H*smoothed_cov*HT;
	TMatrixT<double> invR;
	invertMatrix(R,invR);
	TMatrixT<double> smoothedChiSqu = resT*invR*res;
	return smoothedChiSqu[0][0];
}
