// This Class' Header ------------------
#include "PixHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GeaneTrackRep2.h"
#include "GFDetPlane.h"
#include "GFRectFinitePlane.h"
#include "TRandom.h"
#include "TMath.h"

// Class Member definitions -----------

ClassImp(PixHit)

	TRandom3 PixHit::rand(0);

PixHit::~PixHit()
{}

	PixHit::PixHit()
: PlanarRecoHit(NparHitRep)
{}

PixHit::PixHit(const GFDetPlane& pl,double res) : PlanarRecoHit(NparHitRep){

	fPolicy.setDetPlane(pl);
	fHitCoord[0][0] = rand.Gaus(0.,res);
	fHitCoord[1][0] = rand.Gaus(0.,res);
	fHitCov[0][0] = res*res;
	fHitCov[1][1] = res*res;

}

PixHit::PixHit(TVector3 point,double x_res, double y_res)
	: PlanarRecoHit(NparHitRep){

		fHitCov[0][0] = x_res*x_res;
		fHitCov[1][1] = y_res*y_res;
		GFDetPlane d;

		TVector3 W(0.,0.,1.);
		TVector3 U(1.,0.,0.);
		TVector3 V = W.Cross(U);

//		TVector3 O(point);
		TVector3 O(1.,1.,point.Z());
		d.setO(O);
		d.setU(U);
		d.setV(V);

		point -= O;

		double u,v;
		u = point*U;
		v = point*V;
		fHitCoord[0][0] = u;
		fHitCoord[1][0] = v;
		fPolicy.setDetPlane(d);

}

GFAbsRecoHit* 
PixHit::clone(){
	return new PixHit(*this);
}


	TMatrixT<double>
PixHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
	if(dynamic_cast<const RKTrackRep*>(stateVector) != NULL ) {
		TMatrixT<double> HMatrix(2,5);

		//HMatrix[0][0] = 0.;
		//HMatrix[0][1] = 0.;
		//HMatrix[0][2] = 0.;
		HMatrix(0,3) = 1.;
		//HMatrix[0][4] = 0.;

		//HMatrix[1][0] = 0.;
		//HMatrix[1][1] = 0.;
		//HMatrix[1][2] = 0.;
		//HMatrix[1][3] = 0.;
		HMatrix(1,4) = 1.;
		return HMatrix;
	}
	else if (dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL ) {
    TMatrixT<double> HMatrix(2,6);

    //HMatrix[0][0] = 0.;
    //HMatrix[0][1] = 0.;
    //HMatrix[0][2] = 0.;
    HMatrix(0,3) = 1.;
    //HMatrix[0][4] = 0.;

    //HMatrix[1][0] = 0.;
    //HMatrix[1][1] = 0.;
    //HMatrix[1][2] = 0.;
    //HMatrix[1][3] = 0.;
    HMatrix(1,4) = 1.;
    //HMatrix[1][5] = 0.;
    return HMatrix;
	}
	else {
		std::cerr << "PixHit can only handle state"
			<< " vectors of type GeaneTrackRep2 or RKtrackRep -> abort" 
			<< std::endl;
		throw;
	}
}


