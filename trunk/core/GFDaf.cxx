/* Copyright 2011, Technische Universitaet Muenchen,
   Authors: Karl Bicker, Christian Hoeppner

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

#include "GFDaf.h"
#include "GFTools.h"
#include "GFException.h"
#include "RecoHits/GFDafWireHit.h"
#include "RecoHits/GFAbsWireHit.h"
#include <assert.h>
#include <cmath>


//root stuff
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>
//#define DEBUG


GFDaf::GFDaf(){

	setBetas(81.,8.,4.,1.,1.,1.);
	setProbCut(0.01);
	fKalman.setNumIterations(1);

};

void GFDaf::processTrack(GFTrack* trk) {

#ifdef DEBUG
	std::cout<<"GFDaf::processTrack "<<std::endl;
#endif

	fWeights.clear();

	std::vector<GFDafHit*> eff_hits = initHitsWeights(trk);
	const int nTrkReps = trk->getNumReps();
	const int nDafHits = eff_hits.size();
	if(eff_hits.empty()) {
		for( int i = 0; i != nTrkReps; i++) {
			trk->getTrackRep(i)->setStatusFlag(1);
		}
		return;
	}

	fKalman.initBookkeeping(trk);

	// prepare track
	GFTrack* mini_trk = new GFTrack();
	mini_trk->setSmoothing();

	for( int j=0; j != nDafHits; ++j){
		mini_trk->addHit(eff_hits[j], 0, j); // using dummy det and hit id, they are never used anyway
	}



	// fit for each trackrep separately
	for( int iRep=0; iRep != nTrkReps; ++iRep) { // loop over trackreps

		trk->getBK(iRep)->setNhits(trk->getNumHits());
		if(trk->getTrackRep(iRep)->getStatusFlag()!=0){
			continue;
		}

		mini_trk->addTrackRep(trk->getTrackRep(iRep)); // mini_trk uses the original trackRep here. No copy is made!
		std::vector<std::vector<double> > oldWeights;
		bool oneLastIter = false;
		for( int iBeta=0; iBeta != c_maxIter; ++iBeta) { // loop over. If no convergence is reached after 10 iterations just stop.

#ifdef DEBUG
			std::cout<<"GFDaf::processTrack, trackRep nr. " << iRep << ", beta = " << fBeta[iBeta] << std::endl;
#endif
			for( int j=0; j != nDafHits; j++) {
				GFDafHit* hit = static_cast<GFDafHit*>(mini_trk->getHit(j));
				hit->setWeights(fWeights[iRep][j]);
			}
			if ( iBeta != 0){
				fKalman.blowUpCovs(mini_trk);
			}
			fKalman.processTrack(mini_trk);

			if(mini_trk->getTrackRep(0)->getStatusFlag() != 0){
				break;
			}

			oldWeights = fWeights[iRep];
			try{
				fWeights[iRep] = calcWeights(mini_trk, fBeta[iBeta]);
			} catch(GFException& e) {
				std::cerr<<e.what();
				e.info();
				mini_trk->getTrackRep(0)->setStatusFlag(1);
				break;
			}
			if( oneLastIter == true){
				break;
			}
			if( isConvergent(oldWeights, iRep)){
#ifdef DEBUG
				std::cout << "convergence reached in iteration " << iBeta << "\n";
#endif
				oneLastIter = true;
			}
#ifdef DEBUG
			if( iBeta == c_maxIter-1 ){	std::cout << "giving up after 10 iterations\n";	}
#endif
		} // end loop over betas

		if(mini_trk->getTrackRep(0)->getStatusFlag() == 0 && trk->getSmoothing()){
			copySmoothing(mini_trk, trk, iRep);
		}

		mini_trk->releaseTrackReps();

	} // end loop over trackreps

	saveWeights(trk, mini_trk, fWeights);

	// Calculate the total chi2 and ndf of the DAF fit for all reps with all the stuff that is now in the book keeping
	for(int iRep=0; iRep != nTrkReps; ++iRep) { // loop over trackreps
		GFAbsTrackRep* aTrkRep = trk->getTrackRep(iRep);
		if (aTrkRep->getStatusFlag() != 0){ // not cacluate stuff for track reps that failed during thy
			continue;
		}
		GFBookkeeping* aBK =  trk->getBK(iRep);
		double ndf = 0;
		double totalChi2 = 0;
		double forwardTotalChi2 = 0;

		for( unsigned int i = 0 ; i != trk->getNumHits(); ++i ){
			GFAbsRecoHit* aRecoHit = trk->getHit(i);
			const TMatrixD& H(aRecoHit->getHMatrix(aTrkRep) );
			TVectorD statePredForward = aBK->getVector(GFBKKey_fSt, i);
			TMatrixDSym covPredForward = aBK->getSymMatrix(GFBKKey_fCov, i);
			TVectorD statePredBackward = aBK->getVector(GFBKKey_bSt, i);
			TMatrixDSym covPredBackward = aBK->getSymMatrix(GFBKKey_bCov, i);
			TVectorD w = aBK->getVector(GFBKKey_dafWeight, i);
			TVectorD m(aRecoHit->getRawHitCoord());
			TMatrixDSym V(aRecoHit->getRawHitCov());

			if( m.GetNoElements() == 7){ //we have a cdc Hit.. all this hackery is becaue aRecoHit->getDetPlane crashes sometimes so I cannot use getMeasurment
				double temp = m(6);
				m.ResizeTo(1);
				m(0) = temp;
				temp = V(6,6);
				V.ResizeTo(1,1);
				V(0,0) = temp;
			}
			TMatrixDSym invR(covPredBackward);
			TMatrixDSym invRf(covPredForward);
			invR.Similarity(H);
			invRf.Similarity(H);
			invR += V;
			invRf += V;
			GFTools::invertMatrix(invR);
			GFTools::invertMatrix(invRf);
			if ( w.GetNoElements() == 2) {// this is the special case of when the cdc hit left right ambiguity was resolved by the DAF
				m *= -1; //switch from positive to negative sign
			}
			double unweightedChi2 = invR.Similarity(m - H*statePredBackward);
			double unweightedChi2f = invRf.Similarity(m - H*statePredForward);
			ndf += w(0)*H.GetNrows();
			totalChi2 += w(0)*unweightedChi2;
			forwardTotalChi2 += w(0)*unweightedChi2f;
			if ( w.GetNoElements() == 2) { // this is the special case of when the cdc hit left right ambiguity was resolved by the DAF
				ndf += w(1)*H.GetNrows();
				m *= -1; //switch from negative to positive sign
				unweightedChi2 = invR.Similarity(m - H*statePredBackward);
				unweightedChi2f = invRf.Similarity(m - H*statePredForward);
				totalChi2 += w(1)*unweightedChi2;
				forwardTotalChi2 += w(1)*unweightedChi2f;
			}
		}
		aTrkRep->setNDF(ndf);
		aTrkRep->setChiSqu(totalChi2);
		aTrkRep->setForwardChiSqu(forwardTotalChi2);
	}

	delete mini_trk;


};

std::vector<std::vector<double> > GFDaf::calcWeights(GFTrack* trk, double beta) {

	std::vector<std::vector<double> > ret_val;

	for(unsigned int i=0; i<trk->getNumHits(); ++i) {

		GFDafHit* eff_hit = static_cast<GFDafHit*>(trk->getHit(i));
		unsigned int nEffHits = eff_hit->getNumEffHits();

		std::vector<double> weights;

		if(trk->getBK(0)->hitFailed(i) > 0) { // failed hit
			weights.assign(nEffHits,0.5);
			//std::cout<<"Assumed weight 0.5!!"<<std::endl;
			ret_val.push_back(weights);
			continue;
		}

		std::vector<double> phi;
		double phi_sum = 0;
		double phi_cut = 0;
		TVectorD smoothedState;
		TMatrixDSym smoothedCov;
		GFDetPlane pl;
		GFTools::getBiasedSmoothedData(trk, 0, i, smoothedState, smoothedCov, pl);

		const TMatrixD& H( trk->getHit(i)->getHMatrix(trk->getTrackRep(0)) );
		TVectorD x_smoo(H * smoothedState);



		for(unsigned int j=0; j<nEffHits; j++) {
			double* detV = new double(0);

			try{
				TVectorD m;
				TMatrixDSym Vorig;

				eff_hit->getMeasurement(trk->getTrackRep(0), pl, smoothedState, smoothedCov, m, Vorig, j); // can throw a GFException
				int hitDim = m.GetNoElements();
				TMatrixDSym V( beta * Vorig);
				TVectorD resid(m - x_smoo);
				TMatrixDSym Vinv;
				GFTools::invertMatrix(V, Vinv, detV); // can throw a GFException

				phi.push_back((1./(std::pow(2.*TMath::Pi(),hitDim/2)*sqrt(*detV)))*exp(-0.5*Vinv.Similarity(resid))); // std::pow(double, int) from <cmath> is faster than pow(double, double) from <math.h> when the exponent actually _is_ an integer.
				phi_sum += phi[j];
//				std::cerr << "hitDim " << hitDim << " fchi2Cuts[hitDim] " << fchi2Cuts[hitDim] << std::endl;
				double cutVal = fchi2Cuts[hitDim];
				assert(cutVal>1.E-6);
				//the follwing assumes that in the compeating hits (real hits in one DAF hit) could have different V otherwise calculation could be simplified
				phi_cut += (1./(std::pow(2.*TMath::Pi(),hitDim/2)*sqrt(*detV)))*exp(-0.5*cutVal/beta);
			}
			catch(GFException& e) {
				delete detV;
				e.what();
				e.info();
				phi.push_back(0); //m and Vorig do not contain sensible values, assign weight 0
				continue;
			}

			delete detV;

		}

		for(unsigned int j=0; j<nEffHits; j++) {
			weights.push_back(phi[j]/(phi_sum+phi_cut));
		}

		ret_val.push_back(weights);

	}

	return ret_val;

};

void GFDaf::setProbCut(const double prob_cut){
	for ( int i = 1; i != 6; ++i){
		addProbCut(prob_cut, i);
	}
}

void GFDaf::addProbCut(const double prob_cut, const int measDim){
	if ( prob_cut > 1.0 || prob_cut < 0.0){
		GFException exc("GFDaf::addProbCut prob_cut is not between 0 and 1",__LINE__,__FILE__);
		exc.setFatal();
		throw exc;
	}
	if ( measDim < 1){
		GFException exc("GFDaf::addProbCut measDim must be > 0",__LINE__,__FILE__);
		exc.setFatal();
		throw exc;
	}
	fchi2Cuts[measDim] = ROOT::Math::chisquared_quantile_c( prob_cut, measDim);
}

void GFDaf::setBetas(double b1,double b2,double b3,double b4,double b5,double b6,double b7,double b8,double b9,double b10){
	fBeta.clear();
	assert(b1>0);fBeta.push_back(b1);
	if(b2>0){
		assert(b2<=b1);fBeta.push_back(b2);
		if(b3>=0.) {
			assert(b3<=b2);fBeta.push_back(b3);
			if(b4>=0.) {
				assert(b4<=b3);fBeta.push_back(b4);
				if(b5>=0.) {
					assert(b5<=b4);fBeta.push_back(b5);
					if(b6>=0.) {
						assert(b6<=b5);fBeta.push_back(b6);
						if(b7>=0.) {
							assert(b7<=b6);fBeta.push_back(b7);
							if(b8>=0.) {
								assert(b8<=b7);fBeta.push_back(b8);
								if(b9>=0.) {
									assert(b9<=b8);fBeta.push_back(b9);
									if(b10>=0.) {
										assert(b10<=b9);fBeta.push_back(b10);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	fBeta.resize(c_maxIter,fBeta.back()); //make sure main loop has a maximum of 10 iterations and also make sure the last beta value is used for if more iterations are needed then the ones set by the user.
}

std::vector<GFDafHit*> GFDaf::initHitsWeights(GFTrack* trk) {

	std::vector<GFDafHit*> eff_hits;

	std::vector< std::vector<int> > planes;
	if(not trk->getHitsByPlane(planes)) return eff_hits;
	int nPlanes = planes.size();

	for(int i=0; i<nPlanes; i++) {

		std::vector<GFAbsRecoHit*> hits;

		for(unsigned int j=0; j<planes[i].size(); j++) {
			hits.push_back(trk->getHit(planes[i][j]) );
		}

		GFDafHit* eff_hit;
		if (hits.size()==1 && dynamic_cast<GFAbsWireHit*>(hits[0]) != NULL){
			eff_hit = new GFDafWireHit(dynamic_cast<GFAbsWireHit*>(hits[0]));
		}
		else eff_hit = new GFDafHit(hits);
		eff_hits.push_back(eff_hit);

	}

	for(unsigned int i=0; i<trk->getNumReps(); i++) {
		std::vector<std::vector<double> > rep_weights;
		for(unsigned int j=0; j<eff_hits.size(); j++) {
			rep_weights.push_back(eff_hits[j]->getWeights());
		}
		fWeights.push_back(rep_weights);
	}

	return eff_hits;

}

void GFDaf::copySmoothing(const GFTrack* source, GFTrack* target, int target_irep) {

	unsigned int hit_count = 0;
	GFBookkeeping* sourceBK0 = source->getBK(0);
	GFBookkeeping* targetBKirep = target->getBK(target_irep);

	std::vector< GFBKKey > keys;

	for(unsigned int pl_i=0; pl_i<source->getNumHits(); pl_i++) {

		GFDafHit* eff_hit = static_cast<GFDafHit*>(source->getHit(pl_i));

		for(unsigned int hit_i=0; hit_i<eff_hit->getNumRealHits(); ++hit_i) {

			// copy vectors
			keys = sourceBK0->getVectorKeys();
			for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
				targetBKirep->setVector(keys[iKey], hit_count, sourceBK0->getVector(keys[iKey], pl_i));
			}

			// copy matrices
			keys = sourceBK0->getMatrixKeys();
			for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
				targetBKirep->setMatrix(keys[iKey], hit_count, sourceBK0->getMatrix(keys[iKey], pl_i));
			}

			// copy symmetric matrices
			keys = sourceBK0->getSymMatrixKeys();
			for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
				targetBKirep->setSymMatrix(keys[iKey], hit_count, sourceBK0->getSymMatrix(keys[iKey], pl_i));
			}

			// copy det planes
			keys = sourceBK0->getGFDetPlaneKeys();
			for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
				targetBKirep->setDetPlane(keys[iKey], hit_count, sourceBK0->getDetPlane(keys[iKey], pl_i));
			}

			// copy numbers
			keys = sourceBK0->getNumberKeys();
			for (unsigned int iKey=0; iKey<keys.size(); ++iKey){
				targetBKirep->setNumber(keys[iKey], hit_count, sourceBK0->getNumber(keys[iKey], pl_i));
			}

			++hit_count;
		}
	}

	assert(hit_count == target->getNumHits());
}

void GFDaf::saveWeights(GFTrack* trk, const GFTrack* DafTrack, const std::vector<std::vector<std::vector<double> > >& weights) const {

	assert(weights.size() == trk->getNumReps());

	// get number of real and effective hits
	unsigned int nDafHits(DafTrack->getNumHits());
	std::vector<unsigned int> nRealHits, nEffHits;
	nRealHits.reserve(nDafHits);
	nEffHits.reserve(nDafHits);
	GFDafHit* hit;
	for (unsigned int i=0; i<nDafHits; ++i){
		hit = static_cast<GFDafHit*>(DafTrack->getHit(i));
		nRealHits[i] = hit->getNumRealHits();
		nEffHits[i] = hit->getNumEffHits();
	}

	TVectorD vec;

	for(unsigned int rep_i = 0; rep_i < weights.size(); rep_i++) { // loop over trackReps
		GFBookkeeping* bk = trk->getBK(rep_i);
		bk->bookVectors(GFBKKey_dafWeight);
		unsigned int hit_i = 0;


		for (unsigned int i=0; i<nDafHits; ++i){ // loop over daf hits
			if (nRealHits[i] == nEffHits[i]) {
				for (unsigned int j=0; j<nRealHits[i]; ++j){
					vec.ResizeTo(1);
					vec(0) = weights[rep_i][i][j];
				}
			}
			else {
				assert (nRealHits[i] == 1);
				vec.ResizeTo(nEffHits[i]);
				for (unsigned int j=0; j<nEffHits[i]; ++j){
					vec[j] = weights[rep_i][i][j];
				}
			}
			bk->setVector(GFBKKey_dafWeight, hit_i++, vec);
		}

		assert(hit_i == trk->getNumHits());
	}

}

bool GFDaf::isConvergent(const std::vector<std::vector<double> >& oldWeights, int iRep) const{
	const int n = oldWeights.size();
	const std::vector<std::vector<double> >& newWeights = fWeights[iRep];
	assert(n == int(newWeights.size()));
	for( int i = 0; i != n; ++i){
		const int m = oldWeights[i].size();
		assert(m == int(newWeights[i].size()));
		for( int j = 0; j != m; ++j){
			if( fdim( oldWeights[i][j] , newWeights[i][j]) > 0.001 ){ //Moritz just made the value up. has to be tested if good
				return false;
			}
		}
	}
	return true;
}

