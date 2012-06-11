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
#include<GFDaf.h>

GFDaf::GFDaf() : GFKalman::GFKalman() {

	setBetas(81.,8.,4.,1.,1.,1.);
	setProbCut(0.01);
	GFKalman::setNumIterations(1);

};

void GFDaf::processTrack(GFTrack* trk) {

  fWeights.clear();

	std::vector<GFDafHit*> eff_hits = initHitsWeights(trk);
	if(eff_hits.size() == 0) {
		for(unsigned int i=0; i<trk->getNumReps(); i++) {
			trk->getTrackRep(i)->setStatusFlag(1);
		}
		return;
	}

	GFTrack* mini_trk = new GFTrack();

	for(unsigned int j=0; j<eff_hits.size(); j++) mini_trk->addHit(eff_hits.at(j));

	mini_trk->setSmoothing();

	for(unsigned int i=0; i<trk->getNumReps(); i++) {

		trk->getBK(i)->setNhits(trk->getNumHits());
		if(trk->getTrackRep(i)->getStatusFlag()!=0) continue;

		mini_trk->addTrackRep(trk->getTrackRep(i));

		for(unsigned int iBeta=0; iBeta<fBeta.size(); iBeta++) {

			for(unsigned int j=0; j<mini_trk->getNumHits(); j++) {
				GFDafHit* hit = static_cast<GFDafHit*>(mini_trk->getHit(j));
				hit->setWeights(fWeights.at(i).at(j));
			}			
			if ( iBeta != 0){
				GFKalman::blowUpCovs(mini_trk);
			}
			GFKalman::processTrack(mini_trk);

			if(mini_trk->getTrackRep(0)->getStatusFlag() != 0) break;

			if(iBeta != fBeta.size()-1 ) 
				try{
					fWeights.at(i) = calcWeights(mini_trk, fBeta.at(iBeta));
				} catch(GFException& e) {
					std::cerr<<e.what();
					e.info();
					mini_trk->getTrackRep(0)->setStatusFlag(1);
					break;
				}

		}

		if(trk->getSmoothing()) copySmoothing(mini_trk, trk, i);

		mini_trk->releaseTrackReps();

	}

	saveWeights(trk, fWeights);

	delete mini_trk;

};

std::vector<std::vector<double> > GFDaf::calcWeights(GFTrack* trk, double beta) {

	std::vector<std::vector<double> > ret_val;

	for(unsigned int i=0; i<trk->getNumHits(); i++) {

		GFDafHit* eff_hit = static_cast<GFDafHit*>(trk->getHit(i));
		std::vector<double> phi;
		double phi_sum = 0;
		double phi_cut = 0;

		std::vector<double> weights;

		TMatrixT<double> x_smoo(GFTools::getBiasedSmoothedPos(trk, 0, i));
		TMatrixT<double> smoothedState,smoothedCov;
		GFTools::getBiasedSmoothedData(trk, 0, i, smoothedState,smoothedCov);
		if(x_smoo.GetNrows() == 0) {
			weights.assign(eff_hit->getNumHits(),0.5);
			std::cout<<"Assumed weight 0.5!!"<<std::endl;
			ret_val.push_back(weights);
			continue;
		}

		for(unsigned int j=0; j<eff_hit->getNumHits(); j++) {
			GFAbsRecoHit* real_hit = eff_hit->getHit(j);
			GFDetPlane pl;
			TMatrixT<double> m,Vorig;
			try{
				pl = real_hit->getDetPlane(trk->getTrackRep(0));
				real_hit->getMeasurement(trk->getTrackRep(0),pl,smoothedState,smoothedCov,m,Vorig);
			} catch(GFException& e) {
				std::cerr<<e.what();
				e.info();
				continue; //m and Vorig do not contain sensible values, skip hit
			}

			TMatrixT<double> V( beta * Vorig);
			TMatrixT<double> resid(m - x_smoo);
			TMatrixT<double> resid_T(resid);
			resid_T.T();
			double detV = V.Determinant();
			TMatrixT<double> Vinv;
			GFTools::invertMatrix(V,Vinv);

			phi.push_back((1./(pow(2.*TMath::Pi(),V.GetNrows()/2.)*sqrt(detV)))*exp(-0.5*(resid_T * Vinv * resid)[0][0]));
			phi_sum += phi.at(j);

			double cutVal = fchi2Cuts[V.GetNrows()];
			assert(cutVal>1.E-6);
			phi_cut += (1./(pow(2.*TMath::Pi(),V.GetNrows()/2.)*sqrt(detV)))*exp(-0.5*cutVal/beta);

		}

		for(unsigned int j=0; j<eff_hit->getNumHits(); j++) {

			weights.push_back(phi.at(j)/(phi_sum+phi_cut));

		}
		ret_val.push_back(weights);

	}

	return ret_val;

};

void GFDaf::setProbCut(double val){

	if(fabs(val-0.01)<1.E-10){
		fchi2Cuts[1] = 6.63;
		fchi2Cuts[2] = 9.21;
		fchi2Cuts[3] = 11.34;
		fchi2Cuts[4] = 13.23;
		fchi2Cuts[5] = 15.09;
	}
	else   if(fabs(val-0.005)<1.E-10){
		fchi2Cuts[1] = 7.88;
		fchi2Cuts[2] = 10.60;
		fchi2Cuts[3] = 12.84;
		fchi2Cuts[4] = 14.86;
		fchi2Cuts[5] = 16.75;
	}
	else   if(fabs(val-0.001)<1.E-10){
		fchi2Cuts[1] = 10.83;
		fchi2Cuts[2] = 13.82;
		fchi2Cuts[3] = 16.27;
		fchi2Cuts[4] = 18.47;
		fchi2Cuts[5] = 20.51;
	}
	else{
		GFException exc("GFDafsetProbCut() value is not supported",__LINE__,__FILE__);
		exc.setFatal();
		throw exc;
	}

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
}

std::vector<GFDafHit*> GFDaf::initHitsWeights(GFTrack* trk) {

	std::vector<GFDafHit*> eff_hits;

	std::vector< std::vector<int>* > planes;
	if(not trk->getHitsByPlane(planes)) return eff_hits;
	int nPlanes = planes.size();

	for(int i=0; i<nPlanes; i++) {

		std::vector<GFAbsRecoHit*> hits;

		for(unsigned int j=0; j<planes.at(i)->size(); j++) {
			hits.push_back(trk->getHit(planes.at(i)->at(j)) );
		}

		GFDafHit* eff_hit = new GFDafHit(hits);
		eff_hits.push_back(eff_hit);

	}

	for(unsigned int i=0; i<trk->getNumReps(); i++) {
		std::vector<std::vector<double> > rep_weights;
		for(unsigned int j=0; j<eff_hits.size(); j++) {
			std::vector<double> single_weights;
			single_weights.assign(eff_hits.at(j)->getNumHits(),1.);
			rep_weights.push_back(single_weights);
		}
		fWeights.push_back(rep_weights);
	}

	return eff_hits;

}

void GFDaf::copySmoothing(GFTrack* source, GFTrack* target, int target_irep) {

	for(unsigned int i=0; i<target->getNumReps(); i++) {

		std::vector<std::string> mat_keys = target->getBK(i)->getMatrixKeys();
		bool already_there = false;
		for(unsigned int j=0; j<mat_keys.size(); j++) {
			if(mat_keys.at(j) == "fUpSt") already_there = true;
		}

		if(!already_there) {
			target->getBK(i)->bookMatrices("fUpSt");
			target->getBK(i)->bookMatrices("fUpCov");
			target->getBK(i)->bookMatrices("bUpSt");
			target->getBK(i)->bookMatrices("bUpCov");
			target->getBK(i)->bookGFDetPlanes("fPl");
			target->getBK(i)->bookGFDetPlanes("bPl");
			if(target->getTrackRep(i)->hasAuxInfo()) {
				target->getBK(i)->bookMatrices("fAuxInfo");
				target->getBK(i)->bookMatrices("bAuxInfo");
			}
		}
	}

	int hit_count = 0;

	for(unsigned int pl_i=0; pl_i<source->getNumHits(); pl_i++) {

		GFDafHit* eff_hit = static_cast<GFDafHit*>(source->getHit(pl_i));

		for(unsigned int hit_i=0; hit_i<eff_hit->getNumHits(); hit_i++) {

			TMatrixT<double> fUpSt, fUpCov, bUpSt, bUpCov, fAuxInfo, bAuxInfo;
			GFDetPlane fPl, bPl;

			source->getBK(0)->getMatrix("fUpSt",pl_i,fUpSt);
			source->getBK(0)->getMatrix("fUpCov",pl_i,fUpCov);
			source->getBK(0)->getDetPlane("fPl",pl_i,fPl);
			source->getBK(0)->getMatrix("bUpSt",pl_i,bUpSt);
			source->getBK(0)->getMatrix("bUpCov",pl_i,bUpCov);
			source->getBK(0)->getDetPlane("bPl",pl_i,bPl);

			if(source->getTrackRep(0)->hasAuxInfo()) {
				source->getBK(0)->getMatrix("fAuxInfo",pl_i,fAuxInfo);
				source->getBK(0)->getMatrix("bAuxInfo",pl_i,bAuxInfo);
				target->getBK(target_irep)->setMatrix("fAuxInfo",hit_count,fAuxInfo);
				target->getBK(target_irep)->setMatrix("bAuxInfo",hit_count,bAuxInfo);
			}

			target->getBK(target_irep)->setMatrix("fUpSt",hit_count,fUpSt);
			target->getBK(target_irep)->setMatrix("fUpCov",hit_count,fUpCov);
			target->getBK(target_irep)->setDetPlane("fPl",hit_count,fPl);
			target->getBK(target_irep)->setMatrix("bUpSt",hit_count,bUpSt);
			target->getBK(target_irep)->setMatrix("bUpCov",hit_count,bUpCov);
			target->getBK(target_irep)->setDetPlane("bPl",hit_count,bPl);

			hit_count++;

		}

	}

}

void GFDaf::saveWeights(GFTrack* trk, const std::vector<std::vector<std::vector<double> > >& weights) const {

	assert(weights.size() == trk->getNumReps());

	for(unsigned int rep_i = 0; rep_i < weights.size(); rep_i++) {
		GFBookkeeping* bk = trk->getBK(rep_i);
		bk->bookNumbers("dafWeight");
		unsigned int hit_i = 0;
		for(unsigned int dafhit_i = 0; dafhit_i < weights.at(rep_i).size(); dafhit_i++) {
			for(unsigned int hit_on_pl_i = 0; hit_on_pl_i < weights.at(rep_i).at(dafhit_i).size(); hit_on_pl_i++) {
				bk->setNumber("dafWeight", hit_i, weights.at(rep_i).at(dafhit_i).at(hit_on_pl_i));
				hit_i++;
			}
		}
		assert(hit_i == trk->getNumHits());
	}

}

