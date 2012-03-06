/*
 * G4GFConv.cc
 *
 *  Created on: Feb 2, 2012
 *      Author: poehler
 */

#include "G4GFConv.h"

#include "G4Vector3D.hh"
#include "TVector3.h"
#include "TMatrixT.h"

#include "G4Globals.hh"
#include "G4ErrorSurfaceTrajState.hh"
#include "G4ThreeVector.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorFreeTrajState.hh"
#include <G4ErrorSymMatrix.hh>



TVector3 G4toGFVec (const G4Vector3D& g4vec)
{
	return TVector3(g4vec.x()/10, g4vec.y()/10, g4vec.z()/10);
}

G4Vector3D GFtoG4Vec (const TVector3& gfvec)
{
	return G4Vector3D(gfvec.x()*10, gfvec.y()*10, gfvec.z()*10);
}

void GFcov7toG4cov5 (const TMatrixT<double>& in, G4ErrorTrajErr& out, const TVector3& U, const TVector3& V, const TMatrixT<double>& GFstate7, const double& c)
{
/*	if(in.GetNcols() != 7 || in.GetNrows()() != 7 ){ //|| out.GetNcols() != 5 || out.GetNrows()() != 5 || ){
		std::cerr << "GFcov7toG4cov5 ERROR: Wrong size of matrix\n";
	}*/
	TMatrixT<double> J(5,7);

	double R = GFstate7[6][0];

	double K_u = c*(U.x()*GFstate7[3][0] + U.y()*GFstate7[4][0] + U.z()*GFstate7[4][0]) / R;
	double K_v = c*(V.x()*GFstate7[3][0] + V.y()*GFstate7[4][0] + V.z()*GFstate7[4][0]) / R;



	unsigned int k;

	// 1. col:  0  0  0  0  0  0  c
	for(k = 0; k < 6; k++) J[0][k] = 0;
	J[0][6] = c;

	// 2. col:  0  0  0 c*u_x/R c*u_y/R c*u_z/R  -K_u/R^2

	for(k = 0; k < 3; k++) J[1][k] = 0;
	J[1][3] = c*U.x()/R;
	J[1][4] = c*U.y()/R;
    J[1][5] = c*U.z()/R;
    J[1][6] = (K_u*-1)/(R*R);

	// 3. col:  0  0  0 c*v_x/R c*v_y/R c*v_z/R  -K_v/R^2

	for(k = 0; k < 3; k++) J[2][k] = 0;
	J[2][3] = c*V.x()/R;
	J[2][4] = c*V.y()/R;
    J[2][5] = c*V.z()/R;
    J[2][6] = (K_v*-1)/(R*R);

    // 4. col:  u_x  u_y  u_z  0  0  0  0

   	for(k = 3; k < 6; k++) J[3][k] = 0;
   	J[3][0] = U.X();
   	J[3][1] = U.Y();
   	J[3][2] = U.Z();

    // 5. col:  v_x  v_y  v_z  0  0  0  0

   	for(k = 3; k < 6; k++) J[4][k] = 0;
   	J[4][0] = V.X();
   	J[4][1] = V.Y();
   	J[4][2] = V.Z();




	TMatrixT<double> J_T(7,5);
	J_T = J.T();
	TMatrixT<double> out_(5,5);
	out_ = J*in*J_T;

	for(int i = 0; i < 5; i++){
		for(int j = 0; j < 5; j++)
			out[i][j] = out_[i][j];
	}
}


void G4cov5toGFcov7 (const TMatrixT<double>& in, G4ErrorTrajErr& out, const TVector3& U, const TVector3& V, const TMatrixT<double>& GFstate7, const double& c){

		TMatrixT<double> J(5,7);

		double R = GFstate7[6][0];

		double K_u = c*(U.x()*GFstate7[3][0] + U.y()*GFstate7[4][0] + U.z()*GFstate7[4][0]) / R;
		double K_v = c*(V.x()*GFstate7[3][0] + V.y()*GFstate7[4][0] + V.z()*GFstate7[4][0]) / R;
	//	double K_w =



		unsigned int k;

		// 1. col:  0  0  0  0  0  0  c
		for(k = 0; k < 6; k++) J[0][k] = 0;
		J[0][6] = c;

		// 2. col:  0  0  0 c*u_x/R c*u_y/R c*u_z/R  -K_u/R^2

		for(k = 0; k < 3; k++) J[1][k] = 0;
		J[1][3] = c*U.x()/R;
		J[1][4] = c*U.y()/R;
	    J[1][5] = c*U.z()/R;
	    J[1][6] = (K_u*-1)/(R*R);

		// 3. col:  0  0  0 c*v_x/R c*v_y/R c*v_z/R  -K_v/R^2

		for(k = 0; k < 3; k++) J[2][k] = 0;
		J[2][3] = c*V.x()/R;
		J[2][4] = c*V.y()/R;
	    J[2][5] = c*V.z()/R;
	    J[2][6] = (K_v*-1)/(R*R);

	    // 4. col:  u_x  u_y  u_z  0  0  0  0

	   	for(k = 3; k < 6; k++) J[3][k] = 0;
	   	J[3][0] = U.X();
	   	J[3][1] = U.Y();
	   	J[3][2] = U.Z();

	    // 5. col:  v_x  v_y  v_z  0  0  0  0

	   	for(k = 3; k < 6; k++) J[4][k] = 0;
	   	J[4][0] = V.X();
	   	J[4][1] = V.Y();
	   	J[4][2] = V.Z();


	   	TMatrixT<double> J_inv(7,5);
	   	J_inv = J.Invert();



		TMatrixT<double> J_T(7,5);
		J_T = J.T();
		TMatrixT<double> out_(5,5);
		out_ = J*in*J_T;

		for(int i = 0; i < 5; i++){
			for(int j = 0; j < 5; j++)
				out[i][j] = out_[i][j];
		}



}
