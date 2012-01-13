
#include "GenfitDisplay.h"

#include <assert.h>
#include <cmath>
#include <exception>
#include <iostream>
#include <GFAbsRecoHit.h>
#include <GFAbsTrackRep.h>
#include <GFConstField.h>
#include <GFDetPlane.h>
#include <GFException.h>
#include <GFFieldManager.h>
#include <GFTools.h>
#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveGeoNode.h>
#include <TEveGeoShape.h>
#include <TEveStraightLineSet.h>
#include <TDecompSVD.h>
#include <TGButton.h>
#include <TGeoEltu.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoSphere.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TMatrixDEigen.h>
#include <TROOT.h>
#include <TVector2.h>
#include <TVectorD.h>
#include <TSystem.h>


GenfitDisplay* GenfitDisplay::eventDisplay = NULL;

GenfitDisplay::GenfitDisplay() {

	if((!gApplication) || (gApplication && gApplication->TestBit(TApplication::kDefaultApplication))) {
		std::cout << "In GenfitDisplay ctor: gApplication not found, creating..." << std::flush;
		new TApplication("ROOT_application", 0, 0);
		std::cout << "done!" << std::endl;
	}
	if(!gEve) {
		std::cout << "In GenfitDisplay ctor: gEve not found, creating..." << std::flush;
		TEveManager::Create();
		std::cout << "done!" << std::endl;
	}

	fEventId = 0;
	setOptions();
	setErrScale();

}

void GenfitDisplay::setOptions(std::string opts) { fOption = opts; }

void GenfitDisplay::setErrScale(double errScale) { fErrorScale = errScale; }

double GenfitDisplay::getErrScale() { return fErrorScale; }

GenfitDisplay* GenfitDisplay::getInstance() {

	if(eventDisplay == NULL) {
		eventDisplay = new GenfitDisplay();
	}
	return eventDisplay;

}

GenfitDisplay::~GenfitDisplay() { reset(); }

void GenfitDisplay::reset() { 

	for(unsigned int i = 0; i < fEvents.size(); i++) {

		for(unsigned int j = 0; j < fEvents.at(i)->size(); j++) {

			delete fEvents.at(i)->at(j);
			fEvents.at(i)->clear();
			delete fEvents.at(i);

		}

	}

	fEvents.clear(); 

}

void GenfitDisplay::addEvent(std::vector<GFTrack*>& evts) { 

	std::vector<GFTrack*>* vec = new std::vector<GFTrack*>;

	for(unsigned int i = 0; i < evts.size(); i++) {

		vec->push_back(new GFTrack(*(evts.at(i))));

	}

	fEvents.push_back(vec); 

}

void GenfitDisplay::next(unsigned int stp) { 

	gotoEvent(fEventId + stp);

}

void GenfitDisplay::prev(unsigned int stp) { 

	if(fEventId < (int)stp) {
		gotoEvent(0);
	} else {
		gotoEvent(fEventId - stp);
	}

}

int GenfitDisplay::getNEvents() { return fEvents.size(); }

void GenfitDisplay::gotoEvent(unsigned int id) {

	if(id >= fEvents.size()) id = fEvents.size() - 1;

	fEventId = id;

	std::cout << "At event " << id << std::endl;
	gEve->GetCurrentEvent()->DestroyElements();
	double old_error_scale = fErrorScale;
	drawEvent(fEventId);
	if(old_error_scale != fErrorScale) drawEvent(fEventId); // if autoscaling changed the error, draw again.
	fErrorScale = old_error_scale;

};

void GenfitDisplay::open() {

	bool drawSilent = false;
	bool drawGeometry = false;

	// parse the global options
	for(size_t i = 0; i < fOption.length(); i++) {
		if(fOption.at(i) == 'X') drawSilent = true;
		if(fOption.at(i) == 'G') drawGeometry = true;
	}

	// draw the geometry, does not really work yet. If it's fixed, the docu in the header file should be changed.
	if(drawGeometry) {
		TGeoNode* top_node = gGeoManager->GetTopNode();
		assert(top_node != NULL);

		//Set transparency & color of geometry
		TObjArray* volumes = gGeoManager->GetListOfVolumes();
		for(int i = 0; i < volumes->GetEntries(); i++) {
			TGeoVolume* volume = dynamic_cast<TGeoVolume*>(volumes->At(i));
			assert(volume != NULL);
			volume->SetLineColor(12);
			volume->SetTransparency(50);
		}*/

		TEveGeoTopNode* eve_top_node = new TEveGeoTopNode(gGeoManager, top_node);
		gEve->AddGlobalElement(eve_top_node);
	}

	if(getNEvents() > 0) {
		double old_error_scale = fErrorScale;
		drawEvent(0);
		if(old_error_scale != fErrorScale) gotoEvent(0); // if autoscaling changed the error, draw again.
		fErrorScale = old_error_scale;
	}

	if(!drawSilent) {
		makeGui();
		gApplication->Run(kTRUE);
	}

}

void GenfitDisplay::drawEvent(unsigned int id) {

	// parse the option string ------------------------------------------------------------------------
	bool drawAutoScale = false; 
	bool drawDetectors = false;
	bool drawHits = false;
	bool drawScaleMan = false;
	bool drawTrackMarkers = false;
	bool drawPlanes = false;
	bool drawTrack = false;
	bool drawRawHits = false;

	if(fOption != "") {
		for(size_t i = 0; i < fOption.length(); i++) {
			if(fOption.at(i) == 'A') drawAutoScale = true;
			if(fOption.at(i) == 'D') drawDetectors = true;
			if(fOption.at(i) == 'H') drawHits = true;
			if(fOption.at(i) == 'M') drawTrackMarkers = true;
			if(fOption.at(i) == 'P') drawPlanes = true;
			if(fOption.at(i) == 'S') drawScaleMan = true;
			if(fOption.at(i) == 'T') drawTrack = true;
			if(fOption.at(i) == 'R') drawRawHits = true;
		}
	}
	// finished parsing the option string -------------------------------------------------------------

	// draw SPHits  // quick n dirty hack
	if(drawRawHits){
		for(unsigned int j=0; j<fHits[id].size(); ++j){
			// rotate and translate -------------------------------------------------------
			TGeoMatrix* det_trans = new TGeoGenTrans(fHits[id][j][0], fHits[id][j][1], fHits[id][j][2],
					1./fHits[id][j][3], 1./fHits[id][j][4], 1./fHits[id][j][5], 0);

			TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
			det_shape->SetShape(new TGeoSphere(0.,1.));
			det_shape->SetTransMatrix(*det_trans);
			// finished rotating and translating ------------------------------------------

			det_shape->SetMainColor(kYellow);
			det_shape->SetMainTransparency(70);
			gEve->AddElement(det_shape);
		}
	}





	for(unsigned int i = 0; i < fEvents.at(id)->size(); i++) { // loop over all tracks in an event

		GFTrack* track = fEvents.at(id)->at(i);

		GFAbsTrackRep* rep;
		int irep = 0;
		rep = track->getTrackRep(irep);

		unsigned int numhits = track->getNumHits();
		double charge = rep->getCharge();

		bool smoothing = track->getSmoothing();
		
		if(rep->getStatusFlag()) {
			std::cout << "Warning: Trying to display a track with status flag != 0...";
			if(smoothing) {
				std::cout << "trying without smoothing!";
				smoothing = false;
			}
			std::cout << std::endl;
		}

		TVector3 track_pos;
		TVector3 old_track_pos;

		TEveStraightLineSet* track_lines = NULL;

		// saving the initial state of the representation -----------------------------------------
		GFDetPlane initial_plane = rep->getReferencePlane();
		TMatrixT<double> initial_state(rep->getState());
		TMatrixT<double> initial_cov(rep->getCov());
		TMatrixT<double> initial_auxInfo(*(rep->getAuxInfo(initial_plane)));
		// saved initial state --------------------------------------------------------------------

		for(unsigned int j = 0; j < numhits; j++) { // loop over all hits in the track

			GFAbsRecoHit* hit = track->getHit(j);
			GFDetPlane plane;

			// get the hit infos ------------------------------------------------------------------
			if(smoothing) {
				TMatrixT<double> state;
				TMatrixT<double> cov;
				TMatrixT<double> auxInfo;
				GFTools::getSmoothedData(track, irep, j, state, cov, plane, auxInfo);
				rep->setData(state, plane, &cov, &auxInfo);
			} else {
				try{
					plane = hit->getDetPlane(rep);
					rep->extrapolate(plane);
				}catch(GFException& e) {
					std::cerr << "Error: Exception caught (getDetPlane): Hit " << j << " in Track " << i << " skipped!" << std::endl;
					std::cerr << e.what();
					if (e.isFatal()) {
					  std::cerr<<"Fatal exception, skipping track"<<std::endl;
					  break;
					}
					else continue;
				}
			}
				track_pos = rep->getPos(plane);
			// finished getting the hit infos -----------------------------------------------------

			// sort hit infos into variables ------------------------------------------------------
			TVector3 o = plane.getO();
			TVector3 u = plane.getU();
			TVector3 v = plane.getV();

			std::string hit_type = hit->getPolicyName();

			bool planar_hit = false;
			bool planar_pixel_hit = false;
			bool space_hit = false;
			bool wire_hit = false;
			double_t hit_u = 0;
			double_t hit_v = 0;
			double_t plane_size = 4;
			double_t hit_res_u = 0.5;
			double_t hit_res_v = 0.5;

			TMatrixT<double> hit_coords(hit->getHitCoord(plane));
			TMatrixT<double> hit_cov(hit->getHitCov(plane));
			int hit_coords_dim = hit_coords.GetNrows();

			if(hit_type == "GFPlanarHitPolicy") {
				planar_hit = true;
				if(hit_coords_dim == 1) {
					hit_u = hit_coords(0,0);
					hit_res_u = hit_cov(0,0);
				} else if(hit_coords_dim == 2) {
					planar_pixel_hit = true;
					hit_u = hit_coords(0,0);
					hit_v = hit_coords(1,0);
					hit_res_u = hit_cov(0,0);
					hit_res_v = hit_cov(1,1);
				}
			} else if (hit_type == "GFSpacepointHitPolicy") {
				space_hit = true;
				plane_size = 4; 
			} else if (hit_type == "GFWireHitPolicy") {
				wire_hit = true;
				hit_u = hit_coords(0,0);
				plane_size = 4;
			} else {
				std::cout << "Track " << i << ", Hit " << j << ": Unknown policy name: skipping hit!" << std::endl;
				break;
			}

			if(plane_size < 4) plane_size = 4; 
			// finished setting variables ---------------------------------------------------------

			// draw planes if corresponding option is set -----------------------------------------
			if(drawPlanes || (drawDetectors && planar_hit)) {
				TEveBox* box = boxCreator(track_pos, u, v, plane_size, plane_size, 0.01);
				if (drawDetectors && planar_hit) {
					box->SetMainColor(kCyan);
				} else {
					box->SetMainColor(kGray);
				}
				box->SetMainTransparency(50);
				gEve->AddElement(box);
			}
			// finished drawing planes ------------------------------------------------------------

			// draw track if corresponding option is set ------------------------------------------
			if(drawTrack) {
				if(track_lines == NULL) track_lines = new TEveStraightLineSet;
				if(j > 0) track_lines->AddLine(old_track_pos(0), old_track_pos(1), old_track_pos(2), track_pos(0), track_pos(1), track_pos(2));
				old_track_pos = track_pos;
				if(charge > 0) {
					track_lines->SetLineColor(kRed);
				} else {
					track_lines->SetLineColor(kBlue);
				}
				track_lines->SetLineWidth(2);
				if(drawTrackMarkers) {
					track_lines->AddMarker(track_pos(0), track_pos(1), track_pos(2));
				}
			}
			// finished drawing track -------------------------------------------------------------

			// draw detectors if option is set, only important for wire hits ----------------------
			if(drawDetectors) {

				if(wire_hit) {
					TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
					double pseudo_res_0 = fErrorScale*std::sqrt(hit_cov(0,0));
					if(!drawHits) { // if the hits are also drawn, make the tube smaller to avoid intersecting volumes
						det_shape->SetShape(new TGeoTube(0, hit_u, plane_size));
					} else {
						det_shape->SetShape(new TGeoTube(0, hit_u - pseudo_res_0, plane_size));
					}
					TVector3 norm = u.Cross(v);
					TGeoRotation* det_rot = new TGeoRotation("det_rot",	(u.Theta()*180)/TMath::Pi(), (u.Phi()*180)/TMath::Pi(),
							(norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi(),
							(v.Theta()*180)/TMath::Pi(), (v.Phi()*180)/TMath::Pi()); // move the tube to the right place and rotate it correctly
					TGeoMatrix* det_trans = new TGeoCombiTrans(o(0),o(1),o(2),det_rot);
					det_shape->SetTransMatrix(*det_trans);
					det_shape->SetMainColor(kCyan);
					det_shape->SetMainTransparency(0);
					if((drawHits && (hit_u - pseudo_res_0 > 0)) || !drawHits) {
						gEve->AddElement(det_shape);
					}
				}

			}
			// finished drawing detectors ---------------------------------------------------------

			if(drawHits) {

				// draw planar hits, with distinction between strip and pixel hits ----------------
				if(planar_hit) {
					TVector2 plane_coords = plane.LabToPlane(track_pos);
					if(!planar_pixel_hit) {
						TEveBox* hit_box;
						hit_box = boxCreator((track_pos + (plane_coords.Px() - hit_u)*u), u, v, fErrorScale*std::sqrt(hit_res_u), plane_size, 0.0105);
						hit_box->SetMainColor(kYellow);
						hit_box->SetMainTransparency(0);
						gEve->AddElement(hit_box);
					} else {
						// calculate eigenvalues to draw error-ellipse ----------------------------
						TMatrixDEigen eigen_values(hit_cov);
						TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
						TMatrixT<double> ev = eigen_values.GetEigenValues();
						TMatrixT<double> eVec = eigen_values.GetEigenVectors();
						double pseudo_res_0 = fErrorScale*std::sqrt(ev(0,0));
						double pseudo_res_1 = fErrorScale*std::sqrt(ev(1,1));
						// finished calcluating, got the values -----------------------------------

						// do autoscaling if necessary --------------------------------------------
						if(drawAutoScale) {
							double min_cov = std::min(pseudo_res_0,pseudo_res_1);
							if(min_cov < 1e-5) {
								std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
							} else {
								if(min_cov < 0.049) {
									double cor = 0.05 / min_cov;
									std::cout << "Track " << i << ", Hit " << j << ": Pixel covariance too small, rescaling by " << cor;
									fErrorScale *= cor;
									pseudo_res_0 *= cor;
									pseudo_res_1 *= cor;
									std::cout << " to " << fErrorScale << std::endl; 
								}
							}
						}
						// finished autoscaling ---------------------------------------------------

						// calculate the semiaxis of the error ellipse ----------------------------
						det_shape->SetShape(new TGeoEltu(pseudo_res_0, pseudo_res_1, 0.0105));
						TVector3 pix_pos = track_pos + (plane_coords.Px() - hit_u)*u + (plane_coords.Py() - hit_v)*v;
						TVector3 u_semiaxis = (pix_pos + eVec(0,0)*u + eVec(1,0)*v)-pix_pos;
						TVector3 v_semiaxis = (pix_pos + eVec(0,1)*u + eVec(1,1)*v)-pix_pos;
						TVector3 norm = u.Cross(v);
						// finished calculating ---------------------------------------------------

						// rotate and translate everything correctly ------------------------------
						TGeoRotation* det_rot = new TGeoRotation("det_rot",	(u_semiaxis.Theta()*180)/TMath::Pi(), (u_semiaxis.Phi()*180)/TMath::Pi(),
								(v_semiaxis.Theta()*180)/TMath::Pi(), (v_semiaxis.Phi()*180)/TMath::Pi(),
								(norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi());
						TGeoMatrix* det_trans = new TGeoCombiTrans(pix_pos(0),pix_pos(1),pix_pos(2),det_rot);
						det_shape->SetTransMatrix(*det_trans);
						// finished rotating and translating --------------------------------------

						det_shape->SetMainColor(kYellow);
						det_shape->SetMainTransparency(0);
						gEve->AddElement(det_shape);
					}
				}
				// finished drawing planar hits ---------------------------------------------------

				// draw spacepoint hits -----------------------------------------------------------
				if(space_hit) {

					// get eigenvalues of covariance to know how to draw the ellipsoid ------------
					TMatrixDEigen eigen_values(hit->getRawHitCov());
					TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
					det_shape->SetShape(new TGeoSphere(0.,1.));
					TMatrixT<double> ev = eigen_values.GetEigenValues();
					TMatrixT<double> eVec = eigen_values.GetEigenVectors();
					TVector3 eVec1(eVec(0,0),eVec(1,0),eVec(2,0));
					TVector3 eVec2(eVec(0,1),eVec(1,1),eVec(2,1));
					TVector3 eVec3(eVec(0,2),eVec(1,2),eVec(2,2));
					TVector3 norm = u.Cross(v);
					// got everything we need -----------------------------------------------------


					TGeoRotation* det_rot = new TGeoRotation("det_rot",	(eVec1.Theta()*180)/TMath::Pi(), (eVec1.Phi()*180)/TMath::Pi(),
							(eVec2.Theta()*180)/TMath::Pi(), (eVec2.Phi()*180)/TMath::Pi(),
							(eVec3.Theta()*180)/TMath::Pi(), (eVec3.Phi()*180)/TMath::Pi()); // the rotation is already clear

					// set the scaled eigenvalues -------------------------------------------------
					double pseudo_res_0 = fErrorScale*std::sqrt(ev(0,0));
					double pseudo_res_1 = fErrorScale*std::sqrt(ev(1,1));
					double pseudo_res_2 = fErrorScale*std::sqrt(ev(2,2));
					if(drawScaleMan) { // override again if necessary
						pseudo_res_0 = fErrorScale*0.5;
						pseudo_res_1 = fErrorScale*0.5;
						pseudo_res_2 = fErrorScale*0.5;
					}
					// finished scaling -----------------------------------------------------------

					// autoscale if necessary -----------------------------------------------------
					if(drawAutoScale) {
						double min_cov = std::min(pseudo_res_0,std::min(pseudo_res_1,pseudo_res_2));
						if(min_cov < 1e-5) {
							std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
						} else {
							if(min_cov <= 0.149) {
								double cor = 0.15 / min_cov;
								std::cout << "Track " << i << ", Hit " << j << ": Space hit covariance too small, rescaling by " << cor;
								fErrorScale *= cor;
								pseudo_res_0 *= cor;
								pseudo_res_1 *= cor;
								pseudo_res_2 *= cor;
								std::cout << " to " << fErrorScale << std::endl; 
							}
						}
					}
					// finished autoscaling -------------------------------------------------------

					// rotate and translate -------------------------------------------------------
					TGeoMatrix* det_trans = new TGeoGenTrans(o(0),o(1),o(2),1/(pseudo_res_0),1/(pseudo_res_1),1/(pseudo_res_2),det_rot);
					det_shape->SetTransMatrix(*det_trans);
					// finished rotating and translating ------------------------------------------

					det_shape->SetMainColor(kYellow);
					det_shape->SetMainTransparency(0);
					gEve->AddElement(det_shape);
				}
				// finished drawing spacepoint hits -----------------------------------------------

				// draw wire hits -----------------------------------------------------------------
				if(wire_hit) {
					TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
					double pseudo_res_0 = fErrorScale*std::sqrt(hit_cov(0,0));

					// autoscale if necessary -----------------------------------------------------
					if(drawAutoScale) {
						if(pseudo_res_0 < 1e-5) {
							std::cout << "Track " << i << ", Hit " << j << ": Invalid wire resolution (< 1e-5), autoscaling not possible!" << std::endl;
						} else {
							if(pseudo_res_0 < 0.0049) {
								double cor = 0.005 / pseudo_res_0;
								std::cout << "Track " << i << ", Hit " << j << ": Wire covariance too small, rescaling by " << cor;
								fErrorScale *= cor;
								pseudo_res_0 *= cor;
								std::cout << " to " << fErrorScale << std::endl; 
							}
						}
					}
					// finished autoscaling -------------------------------------------------------

					det_shape->SetShape(new TGeoTube(std::min(0., (double)(hit_u - pseudo_res_0)), hit_u + pseudo_res_0, plane_size));
					TVector3 norm = u.Cross(v);

					// rotate and translate -------------------------------------------------------
					TGeoRotation* det_rot = new TGeoRotation("det_rot",	(u.Theta()*180)/TMath::Pi(), (u.Phi()*180)/TMath::Pi(),
							(norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi(),
							(v.Theta()*180)/TMath::Pi(), (v.Phi()*180)/TMath::Pi());
					TGeoMatrix* det_trans = new TGeoCombiTrans(o(0),o(1),o(2),det_rot);
					det_shape->SetTransMatrix(*det_trans);
					// finished rotating and translating ------------------------------------------

					det_shape->SetMainColor(kYellow);
					det_shape->SetMainTransparency(50);
					gEve->AddElement(det_shape);
				}
				// finished drawing wire hits -----------------------------------------------------

			}

		}

		// reseting to the initial state ----------------------------------------------------------
		rep->setData(initial_state,initial_plane,&initial_cov,&initial_auxInfo);

    try{
      rep->extrapolate(initial_plane);
    }catch(GFException& e) {
      std::cerr << "Error: Exception caught: could not extrapolate back to initial plane " << std::endl;
      std::cerr << e.what();
      continue;
    }
		// done resetting -------------------------------------------------------------------------

		if(track_lines != NULL) gEve->AddElement(track_lines);

	}

	gEve->Redraw3D(kTRUE);

}


void GenfitDisplay::addHits(std::vector<std::vector<double> > hits){
	fHits.push_back(hits);
}




TEveBox* GenfitDisplay::boxCreator(TVector3 o, TVector3 u, TVector3 v, float ud, float vd, float depth) {

	TEveBox* box = new TEveBox;
	float vertices[24];

	TVector3 norm = u.Cross(v);
	u *= (0.5*ud);
	v *= (0.5*vd);
	norm *= (0.5*depth);

	vertices[0] = o(0) - u(0) - v(0) - norm(0);
	vertices[1] = o(1) - u(1) - v(1) - norm(1);
	vertices[2] = o(2) - u(2) - v(2) - norm(2);
	vertices[3] = o(0) + u(0) - v(0) - norm(0);
	vertices[4] = o(1) + u(1) - v(1) - norm(1);
	vertices[5] = o(2) + u(2) - v(2) - norm(2);
	vertices[6] = o(0) + u(0) - v(0) + norm(0);
	vertices[7] = o(1) + u(1) - v(1) + norm(1);
	vertices[8] = o(2) + u(2) - v(2) + norm(2);
	vertices[9] = o(0) - u(0) - v(0) + norm(0);
	vertices[10] = o(1) - u(1) - v(1) + norm(1);
	vertices[11] = o(2) - u(2) - v(2) + norm(2);
	vertices[12] = o(0) - u(0) + v(0) - norm(0);
	vertices[13] = o(1) - u(1) + v(1) - norm(1);
	vertices[14] = o(2) - u(2) + v(2) - norm(2);
	vertices[15] = o(0) + u(0) + v(0) - norm(0);
	vertices[16] = o(1) + u(1) + v(1) - norm(1);
	vertices[17] = o(2) + u(2) + v(2) - norm(2);
	vertices[18] = o(0) + u(0) + v(0) + norm(0);
	vertices[19] = o(1) + u(1) + v(1) + norm(1);
	vertices[20] = o(2) + u(2) + v(2) + norm(2);
	vertices[21] = o(0) - u(0) + v(0) + norm(0);
	vertices[22] = o(1) - u(1) + v(1) + norm(1);
	vertices[23] = o(2) - u(2) + v(2) + norm(2);


	for(int k = 0; k < 24; k += 3) box->SetVertex((k/3), vertices[k], vertices[k+1], vertices[k+2]);

	return box;

}

void GenfitDisplay::makeGui() {

	TEveBrowser* browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kLeft);

	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
	frmMain->SetWindowName("XX GUI");
	frmMain->SetCleanup(kDeepCleanup);

	TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
	{

		TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
		TGPictureButton* b = 0;
		GenfitDisplay*  fh = GenfitDisplay::getInstance(); 

		b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
		hf->AddFrame(b);
		b->Connect("Clicked()", "GenfitDisplay", fh, "prev()");

		b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
		hf->AddFrame(b);
		b->Connect("Clicked()", "GenfitDisplay", fh, "next()");
	}
	frmMain->AddFrame(hf);

	frmMain->MapSubwindows();
	frmMain->Resize();
	frmMain->MapWindow();

	browser->StopEmbedding();
	browser->SetTabTitle("Event Control", 0);
}

ClassImp(GenfitDisplay)
