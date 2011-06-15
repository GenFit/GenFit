/* *************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
$Log: G3Medium.cxx,v $
Revision 1.3  2004/01/28 08:17:52  brun
Reintroduce the Geant3 graphics classes (thanks Andreas Morsch)

Revision 1.1.1.1  2002/07/24 15:56:26  rdm
initial import into CVS

*/


//
// G3 Medium Class for the G3 GUI
// Author: Andreas Morsch
// andreas.morsch@cern.ch
//

#include "G3Medium.h"

ClassImp(G3Medium)

G3Medium::G3Medium()
{
// constructor
    fId=-1;
}

G3Medium::G3Medium(Int_t imed, Int_t imat, const char* name,
			   Int_t isvol, Int_t ifield,
			   Float_t fieldm, Float_t tmaxfd,
			   Float_t stemax, Float_t deemax,
			   Float_t epsil, Float_t stmin)
    : TNamed(name, "Medium")
{
// constructor
    fId=imed;
    fIdMat=imat;
    fIsvol=isvol;
    fIfield=ifield;
    fFieldm=fieldm;
    fTmaxfd=tmaxfd;
    fStemax=stemax;
    fDeemax=deemax;
    fEpsil=epsil;
    fStmin=stmin;
}

Int_t G3Medium::Id()
{
// return medium id
    return fId;
}


Float_t G3Medium::GetPar(Int_t ipar)
{
// Get parameter number ipar
    Float_t p;
    if (ipar < 23) {
	p= fPars[ipar-1];
    } else if(ipar >=23 && ipar <27) {
	p= fPars[ipar-1+3];
    } else {
	p= fPars[ipar-1+4];
    }

    return p;
}

void G3Medium::Streamer(TBuffer &)
{
// dummy streamer
;
}















