/*
 * G4GFConv.hh
 *
 *  Created on: Feb 2, 2012
 *      Author: poehler
 */

#include "G4Vector3D.hh"
#include "TVector3.h"

#include "G4Globals.hh"
#include "G4ErrorSurfaceTrajState.hh"
#include "G4ThreeVector.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorFreeTrajState.hh"
#include <G4ErrorSymMatrix.hh>



#ifndef G4GFCONV_HH_
#define G4GFCONV_HH_

TVector3 G4toGFVec (const G4Vector3D&) ;
G4Vector3D GFtoG4Vec (const TVector3&) ;

void GFcov7toG4cov5 (const TMatrixT<double>&, G4ErrorTrajErr& , const TVector3&, const TVector3&,const TMatrixT<double>&, const double&);
void G4cov5toGFcov7 (const TMatrixT<double>&, G4ErrorTrajErr& , const TVector3&, const TVector3&,const TMatrixT<double>&, const double&);

#endif
