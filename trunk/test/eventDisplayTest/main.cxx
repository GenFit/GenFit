
#include <iostream>
#include <fstream>
#include <string>
#include <TTree.h>
#include <TFile.h>
#include <GFTrack.h>
#include <GFAbsTrackRep.h>
#include <GFDetPlane.h>
#include <GFAbsRecoHit.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TEveBox.h>
#include <TEveManager.h>
#include <TEveBoxSet.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <GFConstField.h>
#include <GFFieldManager.h>
#include <GFException.h>
#include <TEveStraightLineSet.h>
#include <GenfitDisplay.h>

/*#include <TGeoNode.h>
#include <TEveGeoNode.h>
#include <assert.h>*/

int main(int argc, char* argv[]) {


	//TApplication* rootApp = new TApplication("ROOT_application", 0, 0);

	std::string filename = std::string(argv[1]);
	std::string geo_filename = std::string(argv[2]);

	TApplication* rootapp = new TApplication("rootapp", 0, 0);

	TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
	TGeoManager::Import(geo_filename.c_str());
	
	TEveManager::Create();

	GFFieldManager::getInstance()->init(new GFConstField(0.,0.,15.));

	TFile* file = new TFile(filename.c_str());
	TTree* tree = (TTree*)file->Get("o_t");
	TFile* file2 = new TFile("out.root");
	TTree* tree2 = (TTree*)file2->Get("o_t");

	GFTrack* track = new GFTrack();
	GFTrack* track2 = new GFTrack();

	TBranch* branch = tree->GetBranch("GFTrack_branch");
	branch->SetAddress(&track);

	TBranch* branch2 = tree2->GetBranch("GFTrack_branch");
	branch2->SetAddress(&track2);

	//Int_t nentries = tree->GetEntries();
	Int_t nentries = 100;

	GenfitDisplay* display = GenfitDisplay::getInstance();
	display->reset();
	std::vector<GFTrack*> event;

	for(int i = 0; i < nentries; i++) {
		tree->GetEntry(i);
		tree2->GetEntry(i);

		event.clear();

		event.push_back(track);
		event.push_back(track2);

		display->addEvent(event);

	}

	if(argc >= 4) display->setOptions(std::string(argv[3]));

	display->open();
	
	file2->Close();
	file->Close();

}
