#include "SmoothingPulls.h"

void SmoothingPulls::displayPulls(std::vector<std::vector<GFTrack*>*> events, TApplication* rootapp) {

	TCanvas* canvas3 = new TCanvas("canvas3", "Pulls in Plane, Canvas");
	canvas3->Divide(0,2);
	TH1D* histogram7 = new TH1D("Pull_u", "Pull U", 200, -5, 5);
	TH1D* histogram8 = new TH1D("Pull_v", "Pull V", 200, -5, 5);
	TFitResultPtr fit_u;
	TFitResultPtr fit_v;

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptFit(1111);

	std::cout << "--------------------------STARTING GFTOOLS------------------------------" << std::endl;

	for(int i = 0; i < events.size(); i++) {

		for(int j = 0; j < events.at(i)->size(); j++) {
			GFTrack* track = events.at(i)->at(j);
			if(track->getTrackRep(0)->getStatusFlag() != 0 ) {
				std::cout << "Skipping Track " << j << " of event " << i << " because status flag != 0" << std::endl;
				continue;
			}
			for(int k = 0; k < (track->getNumHits()); k++) {
				GFAbsRecoHit* hit = track->getHit(k);
				TMatrixT<double> smoothed_state;
				TMatrixT<double> smoothed_cov;
				GFDetPlane pl; 
				TMatrixT<double> H(hit->getHMatrix(track->getTrackRep(0)));
				try {
					GFTools::getSmoothedData(track, 0, k, smoothed_state, smoothed_cov, pl);
				} catch(GFException& e) {
					std::cerr << e.what();
					std::cerr << "Exception caught in Track " << j << " of event " << i << ". Skipping Track!" << std::endl;
					break;
				}

				TVector2 hit_coords((hit->getHitCoord(pl))(0,0),(hit->getHitCoord(pl))(1,0));
				TVector2 pos((H*smoothed_state)(0,0),(H*smoothed_state)(1,0));
				TVector2 smoothed_res = hit_coords - pos;

				TMatrixT<double> track_cov_tmp(smoothed_cov,TMatrixT<double>::kMultTranspose,H);
			    TMatrixT<double> track_cov(H,TMatrixT<double>::kMult,track_cov_tmp);

				TMatrixT<double> hit_cov(hit->getHitCov(pl));
				TMatrixT<double> smoothed_total_cov(track_cov + hit_cov);

				histogram7->Fill(smoothed_res.X() / (std::sqrt(smoothed_total_cov(0,0))));
				histogram8->Fill(smoothed_res.Y() / (std::sqrt(smoothed_total_cov(1,1))));

			}

		}

	}

	std::cout << "--------------------------STOPPING GFTOOLS------------------------------" << std::endl;

	canvas3->cd(1);
	fit_u = histogram7->Fit("gaus", "S");
	histogram7->Draw();
	canvas3->cd(2);
	fit_v = histogram8->Fit("gaus", "S");
	histogram8->Draw();

}
