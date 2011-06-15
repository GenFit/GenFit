
#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>

#ifndef HELIX_H
#define HELIX_H

class helix {

	public:

	helix();
	helix(TVector3 start, TVector3 dir, TVector3 radius, TVector3 axis);

	void setStartPoint(TVector3 start);
	void setStartDir(TVector3 dir);
	void setRadius(TVector3 Radius);
	void setAxis(TVector3 axis);

	TVector3 getStartPoint() { return fStartPoint; };
	TVector3 getStartDir() { return fStartDir; };
	TVector3 getRadius() { return fRadius; };
	TVector3 getAxis() { return fAxis; };

	bool getChirality(); // true = righthanded
	bool checkConsistency();

	TVector3 getPoint( double t );

	private:

	TVector3 fStartPoint;
	TVector3 fStartDir;
	TVector3 fRadius;
	TVector3 fAxis;

};

#endif
