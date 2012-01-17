
#include <iostream>
#include <GenfitDisplay.h>
#include <GFConstField.h>
#include <GFException.h>
#include <GFFieldManager.h>
#include <GFKalman.h>
#include <GFTools.h>
#include <GFTrack.h>
#include <RKTrackRep.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include "helix.h"
#include "PointHit.h"
#include "SmoothingPulls.h"

int main() {

	TApplication* rootapp = new TApplication("rootapp", 0, 0);
	TEveManager::Create();

	TVector3 old_point(0,0,0);

	double x_res = 0.03;
	double y_res = 0.06;
	double z_res = 0.09;
	double b_field = 10.;

	int err_count = 0;

	TVector3 pos_err(1.,1.,1.);
	TVector3 mom_err(1.,1.,1.);

	TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
	TGeoManager::Import("genfitGeom.root");

	GFFieldManager::getInstance()->init(new GFConstField(0.,0.,b_field));
	TVector3 axis(0,0,1);

	GenfitDisplay* display = GenfitDisplay::getInstance();
	display->reset();

	std::vector<std::vector<GFTrack*>*> events;

	gRandom->SetSeed(4);
//	gRandom->SetSeed(0);

    TCanvas* canvas2 = new TCanvas("canvas2", "Points, Canvas");
    canvas2->Connect("TCanvas", "Closed()", "TApplication", rootapp, "Terminate()");
	canvas2->Divide(0,3);
	TH1D* histogram4 = new TH1D("histogram4", "x", 500, -5, 5);
	TH1D* histogram5 = new TH1D("histogram5", "y", 500, -5, 5);
	TH1D* histogram6 = new TH1D("histogram6", "z", 500, -5, 5);

	for(int i = 0; i < 500; i++) { // looping over events

		std::vector<GFTrack*>* event = new std::vector<GFTrack*>;

		for(int j = 0; j < 1; j++) { // looping over tracks

			TVector3 start(gRandom->Uniform(-10,10),gRandom->Uniform(-10,10),0);
			TVector3 dir(gRandom->Uniform(-10,10),gRandom->Uniform(-10,10),gRandom->Uniform(5,15));
			TVector3 radius = dir.Cross(axis);
			radius.SetMag(gRandom->Uniform(100,200));

			if(gRandom->Uniform(-1,1) > 0) radius *= -1;

			helix hel(start,dir,radius,axis);

			double pos_smear = 0.0001;
			double mom_smear = 0.0001;

			TVector3 start_smeared(gRandom->Gaus(start(0),pos_smear*start(0)),gRandom->Gaus(start(1),pos_smear*start(1)),gRandom->Gaus(start(2),pos_smear*start(2)));

			TVector3 mom_smeared = dir;

			double p_t = (radius.Mag()/100)*(b_field/10)*0.3;
			double dir_t = std::sqrt(std::pow(dir(0),2) + std::pow(dir(1),2));

			mom_smeared.SetMag(dir.Mag()*(p_t/dir_t));

			mom_smeared(0) = gRandom->Gaus(mom_smeared(0),mom_smear*mom_smeared(0));
			mom_smeared(1) = gRandom->Gaus(mom_smeared(1),mom_smear*mom_smeared(1));
			mom_smeared(2) = gRandom->Gaus(mom_smeared(2),mom_smear*mom_smeared(2));

			int pdg_code = 211;
			if (hel.getChirality()) pdg_code *= -1;

			GFAbsTrackRep* rep = new RKTrackRep(start_smeared,mom_smeared,pos_err,mom_err,pdg_code);
			GFTrack* track = new GFTrack(rep, true);

			for(int k = 0; k < 20; k++) { // looping over hits

				TVector3 smeared_point = hel.getPoint(k/80.);
				TVector3 unsmeared_point = smeared_point;

				smeared_point(0) = gRandom->Gaus(smeared_point(0), x_res);
				smeared_point(1) = gRandom->Gaus(smeared_point(1), y_res);
				smeared_point(2) = gRandom->Gaus(smeared_point(2), z_res);

                histogram4->Fill((smeared_point(0) - unsmeared_point(0)) / x_res);
                histogram5->Fill((smeared_point(1) - unsmeared_point(1)) / y_res);
                histogram6->Fill((smeared_point(2) - unsmeared_point(2)) / z_res);

				TVector3 res(x_res,y_res,z_res);
				
				if(smeared_point(0) < 250 && smeared_point(1) < 250 && smeared_point(2) < 250) { 
					track->addHit(new PointHit(smeared_point,res));
				} else {
					std::cout << "Skipped hit " << k << " in track " << j << " of event " << i << ". Out of volume!" << std::endl;
				}

			}

			GFKalman kal;
			std::cout << "At Track " << j << " of Event " << i << "." << std::endl;

			kal.processTrack(track);

			if(track->getTrackRep(0)->getStatusFlag() != 0) err_count++;

			event->push_back(track);

		}

 		display->addEvent(*event);
		events.push_back(event);

	}

	std::cout << "Tracks with errors: " << err_count << std::endl;

	SmoothingPulls::displayPulls(events, rootapp);

	canvas2->cd(1);
	histogram4->Draw();
	canvas2->cd(2);
	histogram5->Draw();
	canvas2->cd(3);
	histogram6->Draw();

	display->setOptions("THDSPM");
	display->open();

}
