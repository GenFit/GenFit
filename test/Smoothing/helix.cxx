
#include "helix.h"

helix::helix() {

	TVector3 str(1,0,0);
	TVector3 dir(0,1,1);
	TVector3 rad(1,0,0);
	TVector3 ax(0,0,1);

	setStartPoint(str);
	setStartDir(dir);
	setRadius(rad);
	setAxis(ax);

}

helix::helix(TVector3 start, TVector3 dir, TVector3 radius, TVector3 axis) {

	setStartPoint(start);
	setStartDir(dir);
	setRadius(radius);
	setAxis(axis);

}

void helix::setStartPoint(TVector3 start) {

	fStartPoint = start;

}

void helix::setStartDir(TVector3 dir) {

	fStartDir = dir;

}

void helix::setRadius(TVector3 radius){

	fRadius = radius;

}

void helix::setAxis(TVector3 axis) {

	axis = axis * (1/axis.Mag()); // normalize axis vector
	fAxis = axis;

}

bool helix::checkConsistency() {

	if((std::abs(fRadius * fStartDir) > 1e-5) || (std::abs(fRadius * fAxis) > 1e-5)) {
		return false;
	} else {
		return true;
	}
}

bool helix::getChirality() {

	if((fRadius.Cross(fStartDir) * fAxis) < 0) {
		return false;
	} else {
		return true;
	}

}

TVector3 helix::getPoint(double t) {


	double t_circ = t;
	double r = fRadius.Mag();
	double a = fStartDir * fAxis;

	t_circ = t*2*TMath::Pi(); 

	double dir = fRadius.Cross(fStartDir) * fAxis;

	if(!checkConsistency()) {
		std::cout << "Parameters not consistent! Aborting..." << std::endl;
		return 0.;
	}

	if(dir > 0) {
		t_circ *= -1; // Check the direction of the rotation 
	} else if(dir < 1e-5 && dir > -1e-5) {
		std::cout << "WARNING: STRAIGHT LINE, DOING NOTHING" << std::endl;
		return 0.;
	}

	TVector3 point(r*std::cos(t_circ)-r, r*std::sin(t_circ), a*t);

	TVector3 r_dir = fRadius * (1/fRadius.Mag());
	TVector3 h = r_dir.Cross(fAxis);

	TMatrixD rotation(3,3);
	rotation(0,0) = r_dir.X();
	rotation(1,0) = r_dir.Y();
	rotation(2,0) = r_dir.Z();
	rotation(0,1) = h.X();
	rotation(1,1) = h.Y();
	rotation(2,1) = h.Z();
	rotation(0,2) = fAxis.X();
	rotation(1,2) = fAxis.Y();
	rotation(2,2) = fAxis.Z();

	point = rotation * point;

	point += fStartPoint;

	return point;

}
