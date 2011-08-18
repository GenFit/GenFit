#include<cmath>

#include"PointHit.h"
#include"PixHit.h"
#include"StripHit.h"

#include<GenfitDisplay.h>

#include<GFAbsTrackRep.h>
#include<GFConstField.h>
#include<GFDaf.h>
#include<GFDetPlane.h>
#include<GFFieldManager.h>
#include<GFKalman.h>
#include<GFTrack.h>

#include<RKTrackRep.h>

#include<TApplication.h>
#include<TCanvas.h>
#include<TGeoManager.h>
#include<TH1D.h>
#include<TRandom.h>
#include<TROOT.h>
#include<TStyle.h>
#include<TVector.h>

int main() {

	TApplication* rootapp = new TApplication("rootapp", 0, 0);

/*	TEveManager::Create();
	GenfitDisplay* display = GenfitDisplay::getInstance();
	display->reset();*/

	int nevs = 500;

	bool strip_hits = true;
	bool pix_hits = true;
	bool noise = true;
	double noise_probability = 0.2;

	double strip_res = 0.005;

	double x_res = 0.03;
	double y_res = 0.06;
	double z_res = 0.09;

	double pos_init_sigma = 0.5;
	double mom_init_sigma = 0.5;

	double dist_btwn_pts = 3;

	double b_field = 10.;

	TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
	TGeoManager::Import("genfitGeom.root");

	GFFieldManager::getInstance()->init(new GFConstField(0.,0.,b_field));

	gRandom->SetSeed(4);

	std::vector<GFTrack*> tracks;
	std::vector<GFTrack*> true_tracks;
	std::vector<std::vector<int> > noise_hits;

	std::vector<TMatrixT<double> > true_states;

	for(int i = 0; i < nevs; i++) { //loop over events (one track events) to generate Tracks

		//Generating Hits with RKTrackRep

		TVector3 starting_pos(0.5,0.,0.);
		TVector3 starting_mom(0.5,1.,1.);
		TVector3 starting_pos_err(1.,1.,1.);
		TVector3 starting_mom_err(1.,1.,1.);

		TVector3 pos_init(gRandom->Gaus(starting_pos.X(),pos_init_sigma),
						  gRandom->Gaus(starting_pos.Y(),pos_init_sigma),
						  gRandom->Gaus(starting_pos.Z(),pos_init_sigma));

		TVector3 mom_init(gRandom->Gaus(starting_mom.X(),mom_init_sigma),
						  gRandom->Gaus(starting_mom.Y(),mom_init_sigma),
						  gRandom->Gaus(starting_mom.Z(),mom_init_sigma));

		GFAbsTrackRep* hitGenerator = new RKTrackRep(starting_pos,starting_mom,starting_pos_err,starting_mom_err,211);

		true_states.push_back(hitGenerator->getState());

		GFAbsTrackRep* rep1 = new RKTrackRep(pos_init, mom_init, starting_pos_err, starting_mom_err, 211);

		pos_init.SetX(gRandom->Gaus(starting_pos.X(),5*pos_init_sigma));
		pos_init.SetY(gRandom->Gaus(starting_pos.Y(),5*pos_init_sigma));
		pos_init.SetZ(gRandom->Gaus(starting_pos.Z(),5*pos_init_sigma));

		mom_init.SetX(gRandom->Gaus(starting_mom.X(),5*mom_init_sigma));
		mom_init.SetY(gRandom->Gaus(starting_mom.Y(),5*mom_init_sigma));
		mom_init.SetZ(gRandom->Gaus(starting_mom.Z(),5*mom_init_sigma));

		GFAbsTrackRep* rep2 = new RKTrackRep(starting_pos, starting_mom, starting_pos_err, starting_mom_err, 211);

		GFTrack* trk = new GFTrack(rep1);
		//trk->addTrackRep(rep2);

		GFTrack* true_trk = new GFTrack(*trk);

		int end = 30;
		int point_end = 10;
		int strip_end = 20;
		int det_id = 0;
		std::vector<int> noise_vector;

		for(int j = 0; j<end; j++) {


			TVector3 pos(hitGenerator->getPos());
			TVector3 mom(hitGenerator->getMom());

			mom.SetMag(dist_btwn_pts);

			GFDetPlane goal_pl(pos+mom,mom);

			TMatrixT<double> state(6,1);
			TMatrixT<double> cov(6,6);

			try{
				hitGenerator->extrapolate(goal_pl,state,cov);
			} catch(GFException& e) {
				e.what();
				std::cerr<<"Error extrapolating in the Hit Generator -> Abort..."<<std::endl;
				exit(1);
			}

			hitGenerator->setData(state, goal_pl, &cov);

			TVector3 hit_pos(hitGenerator->getPos());
			hit_pos.SetX(gRandom->Gaus(hit_pos.X(),x_res));
			hit_pos.SetY(gRandom->Gaus(hit_pos.Y(),y_res));
			if(j < point_end) hit_pos.SetZ(gRandom->Gaus(hit_pos.Z(),z_res));
			TVector3 res(x_res,y_res,z_res);

			if(j < point_end) {
				PointHit* pt_hit = new PointHit(hit_pos, res);
				trk->addHit(pt_hit, 0, j, hit_pos.Mag(), j);
				true_trk->addHit(pt_hit, 0, j, hit_pos.Mag(), j);
				noise_vector.push_back(0);
				if(noise && (gRandom->Uniform(0,1) <= noise_probability)) {
					j++; end++; point_end++; strip_end++;
					hit_pos=hitGenerator->getPos();
					TVector3 noise_pos(gRandom->Uniform(hit_pos.X()-2,hit_pos.X()+2),
							gRandom->Uniform(hit_pos.Y()-2,hit_pos.Y()+2),
							gRandom->Uniform(hit_pos.Z()-2,hit_pos.Z()+2)
							);
					trk->addHit(new PointHit(noise_pos, res), 0, j, hit_pos.Mag(), j);
					noise_vector.push_back(1);
				}
			} else if (j < strip_end && strip_hits) {
				det_id++;
				double strip_res = x_res;
				if(((det_id-1) % 2) == 0) strip_res = y_res;
				StripHit* strp_hit = new StripHit(hit_pos, strip_res, det_id-1);
				trk->addHit(strp_hit, det_id, j, hit_pos.Mag(), det_id);
				true_trk->addHit(strp_hit, det_id, j, hit_pos.Mag(), det_id);
				noise_vector.push_back(0);
				if(noise && (gRandom->Uniform(0,1) <= noise_probability)) {
					j++; end++; strip_end++;
					TVector3 noise_pos(hit_pos);
					if(((det_id-1) % 2) == 0){
						noise_pos.SetY(gRandom->Uniform(hit_pos.Y()-2,hit_pos.Y()+2));
					} else {
						noise_pos.SetX(gRandom->Uniform(hit_pos.X()-2,hit_pos.X()+2));
					}
					trk->addHit(new StripHit(noise_pos, strip_res, det_id-1), det_id, j, hit_pos.Mag(), det_id);
					noise_vector.push_back(1);
				}
			} else if(pix_hits){
				det_id++;
				PixHit* pix_hit = new PixHit(hit_pos, x_res, y_res);
				trk->addHit(pix_hit, det_id, j, hit_pos.Mag(), det_id);
				true_trk->addHit(pix_hit, det_id, j, hit_pos.Mag(), det_id);
				noise_vector.push_back(0);
				if(noise && (gRandom->Uniform(0,1) <= noise_probability)) {
					j++; end++; 
					TVector3 noise_pos(hit_pos);
					noise_pos.SetX(gRandom->Uniform(hit_pos.X()-2,hit_pos.X()+2));
					noise_pos.SetY(gRandom->Uniform(hit_pos.Y()-2,hit_pos.Y()+2));
					trk->addHit(new PixHit(noise_pos, x_res, y_res), det_id, j, hit_pos.Mag(), det_id);
					noise_vector.push_back(1);
				}
			}

		}

		noise_hits.push_back(noise_vector);
		noise_vector.clear();

		tracks.push_back(trk);
		true_tracks.push_back(true_trk);

	}

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);

//ROOT stuff for the histograms
	TCanvas* daf_comp_canvas = new TCanvas("daf_comp_canvas", "Ratio Track Params DAF/True");
	TCanvas* kal_comp_canvas = new TCanvas("kal_comp_canvas", "Ratio Track Params Kalman/True");
	TCanvas* daf_kal_comp_canvas = new TCanvas("daf_kal_comp_canvas", "Ratio Track Params DAF/Kalman");

	daf_comp_canvas->Divide(2,3);
	kal_comp_canvas->Divide(2,3);
	daf_kal_comp_canvas->Divide(2,3);

	TH1D* daf_x1_ratio = new TH1D("daf_x1_ratio", "DAF x[0] - True x[1]", 200, -0.5, 0.5);
	TH1D* daf_x2_ratio = new TH1D("daf_x2_ratio", "DAF x[1] - True x[2]", 200, -0.5, 0.5);
	TH1D* daf_x3_ratio = new TH1D("daf_x3_ratio", "DAF x[2] - True x[3]", 200, -0.5, 0.5);
	TH1D* daf_x4_ratio = new TH1D("daf_x4_ratio", "DAF x[3] - True x[4]", 200, -0.5, 0.5);
	TH1D* daf_x5_ratio = new TH1D("daf_x5_ratio", "DAF x[4] - True x[5]", 200, -0.5, 0.5);
	TH1D* kal_x1_ratio = new TH1D("kal_x1_ratio", "Kalman x[0] - True x[1]", 200, -0.5, 0.5);
	TH1D* kal_x2_ratio = new TH1D("kal_x2_ratio", "Kalman x[1] - True x[2]", 200, -0.5, 0.5);
	TH1D* kal_x3_ratio = new TH1D("kal_x3_ratio", "Kalman x[2] - True x[3]", 200, -0.5, 0.5);
	TH1D* kal_x4_ratio = new TH1D("kal_x4_ratio", "Kalman x[3] - True x[4]", 200, -0.5, 0.5);
	TH1D* kal_x5_ratio = new TH1D("kal_x5_ratio", "Kalman x[4] - True x[5]", 200, -0.5, 0.5);
	TH1D* daf_kal_x1_ratio = new TH1D("daf_kal_x1_ratio", "DAF x[0] - Kalman x[1]", 200, -0.5, 0.5);
	TH1D* daf_kal_x2_ratio = new TH1D("daf_kal_x2_ratio", "DAF x[1] - Kalman x[2]", 200, -0.5, 0.5);
	TH1D* daf_kal_x3_ratio = new TH1D("daf_kal_x3_ratio", "DAF x[2] - Kalman x[3]", 200, -0.5, 0.5);
	TH1D* daf_kal_x4_ratio = new TH1D("daf_kal_x4_ratio", "DAF x[3] - Kalman x[4]", 200, -0.5, 0.5);
	TH1D* daf_kal_x5_ratio = new TH1D("daf_kal_x5_ratio", "DAF x[4] - Kalman x[5]", 200, -0.5, 0.5);

    daf_comp_canvas->Connect("TCanvas", "Closed()", "TApplication", rootapp, "Terminate()");

	//Now for the fitting:

	double missident_noise = 0.;
	double missident_hit = 0.;
	int missident_event = -1;
	double hit_count = 0.;

	for(int i = 0; i<nevs; i++) {

		GFTrack* track = tracks.at(i);
		GFTrack* true_track = true_tracks.at(i);
		GFDaf daf;
		GFKalman kal;
		daf.processTrack(track);
		std::vector<std::vector<std::vector<double> > > weights = daf.getWeights();

		kal.processTrack(true_track);

		for(int j=0; j<track->getNumReps(); j++) {

			int noise_hits_found = 0;
			int real_hits_found = 0;

			hit_count += track->getNumHits();
			int hit_in_ev = 0;

			for(int k=0; k<weights.at(j).size(); k++) {

				for(int l=0; l<weights.at(j).at(k).size(); l++) {
//					std::vector<GFTrack*> event;
					if(weights.at(j).at(k).at(l) < 0.1 && noise_hits.at(i).at(hit_in_ev) == 1) {
						noise_hits_found++;
					} else if(weights.at(j).at(k).at(l) > 0.1 && noise_hits.at(i).at(hit_in_ev) == 0) {
						real_hits_found++;
					} else if(weights.at(j).at(k).at(l) > 0.1 && noise_hits.at(i).at(hit_in_ev) == 1){
//						std::cout<<"NOISE HIT MISSED! Event: "<<i<<", Plane: "<<k<<" Hit: "<<hit_in_ev<<std::endl;
						missident_noise++;
						missident_event = i;
/*						event.push_back(true_track);
						event.push_back(track);*/
					} else if(weights.at(j).at(k).at(l) < 0.1 && noise_hits.at(i).at(hit_in_ev) == 0){
//						std::cout<<"REAL HIT NOISED! Event: "<<i<<", Plane: "<<k<<" Hit: "<<hit_in_ev<<std::endl;
						missident_hit++;
						missident_event = i;
/*				event.push_back(true_track);
				event.push_back(track);*/
					} else {
						std::cout<<"This should not have happened?!"<<std::endl;
					}
					hit_in_ev++;
//					display->addEvent(event);
				}
			}

//			std::cout<<"Noise hits found: " << noise_hits_found << ". Real hits found: " << real_hits_found <<std::endl;

			daf_x1_ratio->Fill( ((track->getTrackRep(j)->getState())[0][0]) - ((true_states.at(i))[0][0]));
			daf_x2_ratio->Fill( ((track->getTrackRep(j)->getState())[1][0]) - ((true_states.at(i))[1][0]));
			daf_x3_ratio->Fill( ((track->getTrackRep(j)->getState())[2][0]) - ((true_states.at(i))[2][0]));
			daf_x4_ratio->Fill( ((track->getTrackRep(j)->getState())[3][0]) - ((true_states.at(i))[3][0]));
			daf_x5_ratio->Fill( ((track->getTrackRep(j)->getState())[4][0]) - ((true_states.at(i))[4][0]));

			kal_x1_ratio->Fill( ((true_track->getTrackRep(j)->getState())[0][0]) - ((true_states.at(i))[0][0]));
			kal_x2_ratio->Fill( ((true_track->getTrackRep(j)->getState())[1][0]) - ((true_states.at(i))[1][0]));
			kal_x3_ratio->Fill( ((true_track->getTrackRep(j)->getState())[2][0]) - ((true_states.at(i))[2][0]));
			kal_x4_ratio->Fill( ((true_track->getTrackRep(j)->getState())[3][0]) - ((true_states.at(i))[3][0]));
			kal_x5_ratio->Fill( ((true_track->getTrackRep(j)->getState())[4][0]) - ((true_states.at(i))[4][0]));

			daf_kal_x1_ratio->Fill( ((track->getTrackRep(j)->getState())[0][0]) - ((true_track->getTrackRep(j)->getState())[0][0]));
			daf_kal_x2_ratio->Fill( ((track->getTrackRep(j)->getState())[1][0]) - ((true_track->getTrackRep(j)->getState())[1][0]));
			daf_kal_x3_ratio->Fill( ((track->getTrackRep(j)->getState())[2][0]) - ((true_track->getTrackRep(j)->getState())[2][0]));
			daf_kal_x4_ratio->Fill( ((track->getTrackRep(j)->getState())[3][0]) - ((true_track->getTrackRep(j)->getState())[3][0]));
			daf_kal_x5_ratio->Fill( ((track->getTrackRep(j)->getState())[4][0]) - ((true_track->getTrackRep(j)->getState())[4][0]));

			for(int m=0;m<5;m++) {

				if((((track->getTrackRep(j)->getState())[m][0]) - ((true_track->getTrackRep(j)->getState())[m][0])) > 0.01) {

					if((std::abs(((track->getTrackRep(j)->getState())[m][0]) - ((true_states.at(i))[m][0]))) >
					   (std::abs(((true_track->getTrackRep(j)->getState())[m][0]) - ((true_states.at(i))[m][0])))) {
						std::cout<<"The winner is DAF!"<<std::flush;
						if(i == missident_event) {
							std::cout<<" And there was a misidentified hit."<<std::endl;
						} else {
							std::cout<<std::endl;
						}
					} else {
						std::cout<<"The winner is Kalman!"<<std::flush;
						if(i == missident_event) {
							std::cout<<" And there was a misidentified hit."<<std::endl;
						} else {
							std::cout<<std::endl;
						}
					}

//					std::cout<<"DAF: "<<((track->getTrackRep(j)->getState())[m][0]) - ((true_states.at(i))[m][0])<<" / Kalman: "<<
//							 ((true_track->getTrackRep(j)->getState())[m][0]) - ((true_states.at(i))[m][0])<<std::endl;

				}

			}
/*				if((((track->getTrackRep(j)->getState())[0][0]) - ((true_track->getTrackRep(j)->getState())[0][0])) > 0.01) {
					std::vector<GFTrack*> event;
					event.push_back(true_track);
					event.push_back(track);
					display->addEvent(event);
				}*/

		}

/*	std::vector<GFTrack*> event;
//	event.push_back(true_track);
event.push_back(tracks.at(i));
//	event.push_back(track);
	display->addEvent(event);*/

	}

	std::cout<<"Wrongly identified hits: "<<(missident_noise+missident_hit)<<" of "<<hit_count<<" ("<<(((missident_noise+missident_hit)/hit_count)*100.)<<"%)"<<std::endl;
	std::cout<<"Noise wrongly considered as hit: "<<missident_noise<<std::endl;
	std::cout<<"Hit wrongly considered as noise: "<<missident_hit<<std::endl;

	daf_comp_canvas->cd(1);
	daf_x1_ratio->Draw();
	daf_comp_canvas->cd(2);
	daf_x2_ratio->Draw();
	daf_comp_canvas->cd(3);
	daf_x3_ratio->Draw();
	daf_comp_canvas->cd(4);
	daf_x4_ratio->Draw();
	daf_comp_canvas->cd(5);
	daf_x5_ratio->Draw();

	kal_comp_canvas->cd(1);
	kal_x1_ratio->Draw();
	kal_comp_canvas->cd(2);
	kal_x2_ratio->Draw();
	kal_comp_canvas->cd(3);
	kal_x3_ratio->Draw();
	kal_comp_canvas->cd(4);
	kal_x4_ratio->Draw();
	kal_comp_canvas->cd(5);
	kal_x5_ratio->Draw();

	daf_kal_comp_canvas->cd(1);
	daf_kal_x1_ratio->Draw();
	daf_kal_comp_canvas->cd(2);
	daf_kal_x2_ratio->Draw();
	daf_kal_comp_canvas->cd(3);
	daf_kal_x3_ratio->Draw();
	daf_kal_comp_canvas->cd(4);
	daf_kal_x4_ratio->Draw();
	daf_kal_comp_canvas->cd(5);
	daf_kal_x5_ratio->Draw();

/*	display->setOptions("THDSPM");
	display->open();*/

	rootapp->Run();

}
