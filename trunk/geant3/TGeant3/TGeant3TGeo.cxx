
/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
$Log: TGeant3TGeo.cxx,v $
Revision 1.17  2006/12/19 13:16:19  brun
from Mohammad Al-Turany & Denis Bertini

Changes in  TGeant3/TGeant3.cxx and TGeant3.h
------------------------------------
1. Geane interface functions are added:
    void  eufill(Int_t n, Float_t *ein, Float_t *xlf);
    void  eufilp(const int n,Float_t *ein,Float_t *pli,Float_t *plf);
    void  eufilv(Int_t n, Float_t *ein,	Char_t *namv, Int_t *numv,Int_t *iovl);
    void  trscsd(Float_t *pc,Float_t *rc,Float_t *pd,Float_t *rd,
		 Float_t *h,Float_t *ch,Int_t *ierr,Float_t *spu,Float_t *dj,Float_t *dk);
    void  trsdsc(Float_t *pd,Float_t *rd,Float_t *pc,Float_t *rc,
                           Float_t *h,Float_t *ch,Int_t *ierr,Float_t *spu,Float_t *dj,Float_t *dk);
    void  trscsp(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,
  			Float_t *ch,Int_t *ierr,Float_t *spx);
    void  trspsc(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,
  			Float_t *ch,Int_t *ierr,Float_t *spx);

2. The Gfang function wrapper is added
    void  g3fang( Float_t *, Float_t &,Float_t &, Float_t &, Float_t &,Int_t & );




changes in TGeant3/TGeant3gu.cxx
--------------------------
Adding GCalor interface
	1. function calsig() and gcalor() are added
	2. setting ihadr=5 will call the GCalor routine




changes in TGeant3/TGeant3.h
-----------------------
1. Structures for Geane output are setted as public so that the user can access them

    Ertrio_t *fErtrio
    Eropts_t *fEropts
    Eroptc_t *fEroptc
    Erwork_t *fErwork
    Trcom3_t *fTrcom3

2. The size of the error matrix errin is corrected to 15

Revision 1.15  2006/05/30 13:39:07  brun
From Andrei Gheata:
a patch cleaning-up the usage of TGeant3::Vname method inside TGeant3TGeo.cxx. The truncation of names is not needed when working with TGeo. In most of the cases the names were truncated to 4 chars but the result was not used in the subsequent call to TGeoMCGeometry (this is why it still worked) but a cleanup is good anyway...

Revision 1.14  2006/05/23 15:53:11  brun
From Ivana:
 Adding CurrentVolPath() overloading TGeant3 implementation
 which does not work correctly with longer volume names
 (Oleg Yushchenko)

Revision 1.13  2005/11/18 21:25:22  brun
From Bjorn, Andrei:
Implemented new VMC functions for access to geometry;
added -Woverloaded-virtual to Makefile.linux

Revision 1.12  2005/07/28 12:02:17  brun
From Andrei:
- Fixed problem of material indices when coming from FORTRAN code
(TGeant3TGeo::Gsmate)
- when loading the geometry from file, the value of radlen stored in TGeo is
injected in G3 (with negative sign not to be recomputed by G3) but abslen is
recomputed since it is not valid in TGeo.

Revision 1.11  2005/07/27 13:06:52  brun
Simplify logic in TGeant3TGeo::Mixture (Float_t* case)

Revision 1.10  2005/07/21 17:54:37  brun
Implement same code in the float* versions that were already implemented
in the Double* versions.

Revision 1.9  2005/07/20 09:22:51  brun
From Federico:
Fixes to compile with gcc4CVS: ----------------------------------------------------------------------

Revision 1.8  2005/07/13 09:36:18  brun
From Federico:
  Mods for Mac and removal of stupid printout.

Revision 1.7  2005/07/11 12:00:47  brun
From Andrei Gheata:
a fix in TGeant3TGeo::Gsmixt. The problem (found by Federico) was
that the array of weights for mixture components in case nlmat<0 (number of
atoms) was recomputed once by G3 itself then again inside Gsmixt, resulting
in wrong fractions.

Revision 1.6  2005/06/15 08:49:21  brun
From Andrei Gheata:
Change related to  usage of assemblies (that have no medium).

Revision 1.5  2005/05/17 12:48:00  brun
From Ivana:
- Set name TGeant3TGeo

Revision 1.4  2005/05/11 11:46:54  brun
From Andrei Gheata:
 Implementation of gtmany (PUSH) and glvolu (POP) in
TGeant3TGeo class (needed for the correct boundary crossing sequence in
case of Cerenkov transport).

Revision 1.3  2005/02/08 11:22:03  brun
From Ivana:
For TGeant3.h:
Added IsRootGeometrySupported() function
(now required by TVirtualMC)

For TGeant3.cxx:
Updated text in Fatal in SetRootGeometry.

Revision 1.2  2004/12/21 15:34:48  brun
Implement TGeant3TGeo::isRootGeometry returning kTRUE

Revision 1.1  2004/12/17 11:55:47  brun
A new class TGeant3TGeo (deriving from TGeant3) is introduced.
TGeant3 uses by default the geant3 geometry. TGeant3TGeo uses the TGeo classes.
The two classes are built in the same library. The choice of which version to use
can now be made dynamically at run time (eg based on an environment variable)
or a job control option.
For example the examples like gexam1, gexam4 have been modified to run
with either geant3 geometry or TGeo. To run gexam1 with geant3 geometry do
  gexam1 TGeant3
to run with TGeo (default) do
  gexam1
The examples E01.C, E02 and E03 have also been modified to select the option
at run time based on the environment variable TVirtualMC (set in .rootrc)
If TVirtualMC is set to TGeant3TGeo TGeo geometry will be used.
The Makefile has been modified to build the two classes in the same library.

Revision 1.37  2004/11/23 14:52:52  brun
From Andreas Morsch:
on request of ALICE/TPC I added a new method to TVirtualMC.h

void TVirtualMC::ForceDecayTime(Float_t);


This allows to force the decay time of the current particle.
The use-case implemented in AliRoot is the decay of primary particles
within a user defined radius range.

Revision 1.36  2004/10/13 10:38:32  brun
From Andrei Gheata:
some modifications in the current TGeant3.cxx :

from Mihaela:
- modifications in STATISTICS option: added branches to statistics tree
(statsame + statpath)
- global gckine added. gckine->itrtyp == 7 used in gtnext to optimize
speed - 1% gain (computation of global matrix only when called from gtckov)
- bias of 1E-7 (used previously for making sure a boundary is crossed)
eliminated

I removed the option WITHBOTH and fixed a problem in VolId() - in the
last version of AliRoot some detector was calling gMC->VolId(name) with
a name containing a blank at the end and now all volumes have the blanks
supressed. Fixing this I noticed that there are 4 detectors that in
their StepManager() they call at each step things like:
  if (gMC->CurrentVolId(copy) == gMC->VolId("RICH")) ...
Incredible !!! In G3 native this search by name does not penalize so
much since names are converted to Int_t and the volume bank is looked
for this Int_t. In TGeo we cannot do this since we support long names so
we have to go to gGeoManager->GetVolume("name") which scans a list of
2000 objects at each step several times...   I fixed this by hand by
puting static variables in these methods and Peter will commit the
changes. The gain in speed for TGeo case is considerable with full AliRoot.

Revision 1.35  2004/10/12 07:46:23  brun
>From Ivana:
Implemented new functions from TVirtualMC:
  Int_t NofVolDaughters(const char* volName) const;
  const char*  VolDaughterName(const char* volName, Int_t i) const;
  Int_t        VolDaughterCopyNo(const char* volName, Int_t i) const;
  const char* CurrentVolPath();

Revision 1.34  2004/09/17 08:51:55  brun
>From Ivana
 SetRootGeometry() allowed only with WITHROOT option;
 added Fatal() for other modes.

Revision 1.33  2004/08/25 07:28:54  brun
>From Ivana and Lionel Chaussard
 In method DefineParticles(), we had:

         pdgcode(Tau+)=15
         pdgcode(Tau-)=-15

 This is in contradiction with the other leptons (pdgcode>0 for
 negative leptons, pdgcode<0 for positive leptons). I also checked
 in the PDG WEB pages that Tau- should have a code +15.

 changed to:
         pdgcode(tau+)=-15
         pdgcode(tau-)=+15 ?

Revision 1.32  2004/08/05 12:20:39  brun
>From Andrei Gheata:
I have found/fixed a bug in TGeoManager::IsSameLocation(x,y,z). Also I
have eliminated the penalizing check of IsSameLocation() in gtnext().

Revision 1.31  2004/07/09 12:15:12  brun
>From Ivana:
in case a user defines geometry via TGeo and associates more tracking
 media with the same material, TGeant3 duplicates this material
 for each tracking medium.
 I haven't found a function for getting the number of
 materials/media (TList) so I count them by a loop
 - maybe it can be done more intelligently...

Revision 1.30  2004/07/09 08:11:29  brun
Fix by Ivana/Andrei to call TGeoMedium::setId and not TGeoMedium::SetUniqueID

Revision 1.29  2004/06/17 13:56:53  rdm
changed several "const int" arguments to "int". Was causing warnings of
type "qualifier is meaningless".

Revision 1.28  2004/06/08 10:27:19  brun
>From Ivana:
- Added Bool_t return value to methods
  SetCut(), SetProcess(), DefineParticle(), DefineIon()
- Removed  DefineParticles()

Revision 1.27  2004/05/28 13:45:00  brun
>From Ivana
Implementation of StopRun (new function in TVirtualMC)

Revision 1.26  2004/05/14 08:32:01  brun
In function gtnext, call GetNextBoundary(-step) instead of (step).
This fixes a problem when tracking Cherenkov photons.
(Thanks to Yuri Kharlov for reporting the problem and Andrei for fixing it)

Revision 1.25  2004/03/23 11:16:44  brun
>From Ivana
With the previous changes by Andrei, all fixes by Ivana were lost.
This patch merges Ivana and Andrei versions.

Revision 1.23  2004/03/15 12:18:45  brun
>From Andrei Gheata:
 - minor modifications to cope with geometry retreival from file:
 - ConstructGeometry does not need to be called
 - CloseGeometry not needed

Revision 1.22  2004/02/03 12:47:34  brun
>From Andrei Gheata:
TGeant3:

- calls to gtonly return now always a true value (G3 is seeing an ONLY
geometry with TGeo)
- IsSameLocation inside gtnext not yet eliminated, but I am getting only
1 exception instead of 400 now (when the location really changes)

Revision 1.21  2004/01/28 18:05:24  brun
New version from Peter Hristov adding the graphics interface

Revision 1.20  2004/01/28 08:30:54  brun
Change the call to TRandom::RndmArray in function grndm

Revision 1.19  2004/01/28 08:14:48  brun
Add a CPP option STATISTICS to monitor the fequency of calls to the geometry functions.

Revision 1.18  2003/12/10 15:39:37  brun
iFollowing recent improvements by Andrei, replace in ggperp
the computation of normals:
    Double_t *dblnorm = gGeoManager->FindNormal(kFALSE);
with :
    Double_t *dblnorm = gGeoManager->FindNormalFast();

Revision 1.17  2003/12/10 10:32:09  brun
Add a protection in TGeant3TGeo::Gsmate in case the material density is null

Revision 1.16  2003/12/01 23:51:22  brun
>From Andrei and Peter:
add a few missing cases when compiling with the WITHROOT option.

Revision 1.15  2003/11/28 09:44:15  brun
New version of TGeant3 supporting the options WITHG3 and WITHROOT

Revision 1.14  2003/10/09 06:28:45  brun
In TGeant3TGeo::ParticleName, increase size of local array name[20] to name[21]

Revision 1.13  2003/09/26 15:01:08  brun
>From Ivana;
- implemented new functions from TVirtualMC
  enabling user to define own particles and ions
  + getter functions::
    DefineParticle(..)
    DefineIon(..)
    ParticleName(..) const
    ParticleMass(..) const
    Double_t  ParticleCharge(..) const
    Double_t  ParticleLifeTime(..) const
- corrected charge in AddParticlesToPdgDataBase

Revision 1.12  2003/07/22 06:53:28  brun
This version does not yet support TGeo geometry.
TVirtualMC must be initialized with the 3rd argument set to kFALSE

Revision 1.11  2003/07/18 10:22:50  brun
Changes to reflect the equivalent changes in the abstract classes in vmc
(thanks Peter Hristov)

Revision 1.10  2003/07/16 07:40:09  brun
>From Andreas Morsch

- default g3 specific initialisation moved to TGeant3TGeo::Init()
  (This avoids the cast to TGeant3* in the Config.C)
- "CKOV" added to SetProcess

Revision 1.9  2003/06/03 21:26:46  brun
New version of gustep by Andreas Morsch

Revision 1.8  2003/02/28 10:41:49  brun
>From Andreas
 In DefineParticles(): rho0 decay channel corrected

Revision 1.7  2003/02/04 17:50:34  brun
>From Ivana
 In Mixture(): pass abs(nlmat) to CreateFloatArray calls
 as nlmat can be negative.

Revision 1.6  2003/01/31 18:23:06  brun
Ivana suggested corrections.
- corrected tau pdg code
- Warning if external decayer needed but not defined.

Revision 1.5  2003/01/23 11:34:04  brun
In gustep, replace
   gMC->TrackPosition(x,y,z);
by
   geant3->TrackPosition(x,y,z);

Revision 1.4  2003/01/06 17:20:52  brun
Add new functions TrackPosition and TrackMomentum as alternative to the original
functions filling a TLorentzVector object.
Use these new functions in gustep and gudcay.
This makes a 25 per cent speed improvement in case of Alice.

Revision 1.3  2002/12/10 07:58:36  brun
Update by Federico for the calls to Grndm

Revision 1.2  2002/12/06 16:50:30  brun
>From Federico:
the following modifications provide an >6% improvement in speed for
AliRoot.

Revision 1.1.1.1  2002/07/24 15:56:26  rdm
initial import into CVS

Revision 1.5  2002/07/10 09:33:19  hristov
Array with variable size created by new

Revision 1.4  2002/07/10 08:38:54  alibrary
Cleanup of code

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Interface Class to the Geant3.21 MonteCarlo                              //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <ctype.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TArrayI.h"
#include "TArrayD.h"

#include "TGeant3TGeo.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMCGeometry.h"

#include "TCallf77.h"
#include "TVirtualMCDecayer.h"
#include "TPDGCode.h"

#ifndef WIN32
# define g3smate  g3smate_
# define g3smixt  g3smixt_
# define g3stmed  g3stmed_
# define g3treve  g3treve_
# define gtreveroot  gtreveroot_
# define gcomad gcomad_

# define g3brelm g3brelm_
# define g3prelm g3prelm_

#else

# define gzebra  GZEBRA
# define grfile  GRFILE
# define gpcxyz  GPCXYZ
# define ggclos  GGCLOS
# define ginit   GINIT
# define g3cinit  G3CINIT
# define grun    GRUN
# define gtrig   GTRIG
# define gtrigc  GTRIGC
# define gtrigi  GTRIGI
# define gwork   GWORK
# define gfmate  GFMATE
# define gfpart  GFPART
# define gftmed  GFTMED
# define gftmat  GFTMAT
# define gsmate  GSMATE
# define gsmixt  GSMIXT
# define gstmed  GSTMED
# define gstpar  GSTPAR
# define gfkine  GFKINE
# define gfvert  GFVERT
# define gskine  GSKINE
# define gsvert  GSVERT
# define gphysi  GPHYSI
# define gekbin  GEKBIN
# define gfinds  GFINDS
# define gsking  GSKING
# define gskpho  GSKPHO
# define gsstak  GSSTAK
# define gtreve  GTREVE
# define gtreveroot  GTREVEROOT
# define gdtom   GDTOM
# define gmedia  GMEDIA
# define gmtod   GMTOD
# define gsdvn   GSDVN
# define gsdvn2  GSDVN2
# define gsdvs   GSDVS
# define gsdvs2  GSDVS2
# define gsdvt   GSDVT
# define gsdvt2  GSDVT2
# define gsord   GSORD
# define gspos   GSPOS
# define gsposp  GSPOSP
# define gsrotm  GSROTM
# define gprotm  GPROTM
# define gsvolu  GSVOLU
# define gprint  GPRINT
# define gcomad  GCOMAD

#endif

//____________________________________________________________________________
extern "C"
{
  //
  // Prototypes for GEANT functions
  //

  void type_of_call g3treve();

  void type_of_call gtreveroot();

  void type_of_call g3smate(const Int_t&, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t &, Float_t *,
			   Int_t & DEFCHARL);

  void type_of_call g3smixt(const Int_t&, DEFCHARD, const Float_t *,
               const Float_t *, const Float_t &, const Int_t &,
               Float_t * DEFCHARL);

  void type_of_call g3stmed(const Int_t&, DEFCHARD, Int_t &, Int_t &, Int_t &,
			   Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL);

  void type_of_call gcomad(DEFCHARD, Int_t*& DEFCHARL);
}


#ifndef WIN32
#  define gtreveroot gtreveroot_

#else
#  define gtreveroot GTREVEROOT

#endif

extern "C" type_of_call void gtreveroot();
extern "C" type_of_call void ginvolTGeo(Float_t*, Int_t&);
extern "C" type_of_call void gtmediTGeo(Float_t*, Int_t&);
extern "C" type_of_call void gtmanyTGeo(Int_t&);
extern "C" type_of_call void gtonlyTGeo(Int_t&);
extern "C" type_of_call void gmediaTGeo(Float_t*, Int_t&, Int_t&);
extern "C" type_of_call void glvoluTGeo(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier);
extern "C" type_of_call void gtnextTGeo();
extern "C" type_of_call void ggperpTGeo(Float_t*, Float_t*, Int_t&);

//
// Geant3 global pointer
//
Gcvol1_t *gcvol1 = 0;
TGeoNode *gCurrentNode = 0;
R__EXTERN Gctrak_t *gctrak;
R__EXTERN Gcvolu_t *gcvolu;
R__EXTERN Gckine_t *gckine;
R__EXTERN TGeant3* geant3;
R__EXTERN Gcchan_t *gcchan;
R__EXTERN Int_t count_ginvol;
R__EXTERN Int_t count_gmedia;
R__EXTERN Int_t count_gtmedi;
R__EXTERN Int_t count_gtnext;

R__EXTERN void (*fginvol)(Float_t*, Int_t&);
R__EXTERN void (*fgtmedi)(Float_t*, Int_t&);
R__EXTERN void (*fgtmany)(Int_t&);
R__EXTERN void (*fgtonly)(Int_t&);
R__EXTERN void (*fgmedia)(Float_t*, Int_t&, Int_t&);
R__EXTERN void (*fglvolu)(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier);
R__EXTERN void (*fgtnext)();
R__EXTERN void (*fggperp)(Float_t*, Float_t*, Int_t&);


//____________________________________________________________________________
TGeant3TGeo::TGeant3TGeo()
  : TGeant3(),
    fMCGeo(0)
{
  //
  // Default constructor
}

//____________________________________________________________________________
TGeant3TGeo::TGeant3TGeo(const char *title, Int_t nwgeant)
       : TGeant3(title,nwgeant),
         fImportRootGeometry(kFALSE)
{
  //
  // Standard constructor for TGeant3 with ZEBRA initialisation
  //

  SetName("TGeant3TGeo");

  fMCGeo = new TGeoMCGeometry("MCGeo", "TGeo Implementation of VirtualMCGeometry");

  LoadAddress();
  //set pointers to tracker functions
  fginvol = ginvolTGeo;
  fgtmedi = gtmediTGeo;
  fgtmany = gtmanyTGeo;
  fgtonly = gtonlyTGeo;
  fgmedia = gmediaTGeo;
  fglvolu = glvoluTGeo;
  fgtnext = gtnextTGeo;
  fggperp = ggperpTGeo;
}

//____________________________________________________________________________
TGeant3TGeo::~TGeant3TGeo()
{
   delete fMCGeo;
}

//____________________________________________________________________________
void TGeant3TGeo::LoadAddress()
{
  //
  // Assigns the address of the GEANT common blocks to the structures
  // that allow their access from C++
  //
//   printf("LoadAddress\n");
//   TGeant3::LoadAddress();
   gcomad(PASSCHARD("GCVOL1"),(int*&) fGcvol1  PASSCHARL("GCVOL1"));
   gcvol1 = fGcvol1;
}

//_____________________________________________________________________________
void TGeant3TGeo::GeomIter()
{
  //
  // Geometry iterator for moving upward in the geometry tree
  // Initialise the iterator
  //
  fNextVol=gGeoManager->GetLevel();
}

//____________________________________________________________________________
Int_t TGeant3TGeo::NextVolUp(Text_t *name, Int_t &copy)
{
  //
  // Geometry iterator for moving upward in the geometry tree
  // Return next volume up
  //
  fNextVol--;
  if (fNextVol>=0) {
     Int_t level = gGeoManager->GetLevel();
     if (level<=fNextVol) return 0;
     TGeoNode *mother = gGeoManager->GetMother(level-fNextVol);
     if (!mother) return 0;
     sprintf(name, "%s", mother->GetVolume()->GetName());
     copy = mother->GetNumber();
     return mother->GetVolume()->GetNumber();
  }
  return 0;
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::CurrentVolID(Int_t &copy) const
{
  //
  // Returns the current volume ID and copy number
  //
  if (gGeoManager->IsOutside()) return 0;
  TGeoNode *node = gGeoManager->GetCurrentNode();
  copy = node->GetNumber();
  Int_t id = node->GetVolume()->GetNumber();
  return id;
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::CurrentVolOffID(Int_t off, Int_t &copy) const
{
  //
  // Return the current volume "off" upward in the geometrical tree
  // ID and copy number
  //
  if (off<0 || off>gGeoManager->GetLevel()) return 0;
  if (off==0) return CurrentVolID(copy);
  TGeoNode *node = gGeoManager->GetMother(off);
  if (!node) return 0;
  copy = node->GetNumber();
  return node->GetVolume()->GetNumber();
}

//_____________________________________________________________________________
const char* TGeant3TGeo::CurrentVolName() const
{
  //
  // Returns the current volume name
  //
  if (gGeoManager->IsOutside()) return gGeoManager->GetTopVolume()->GetName();
  return gGeoManager->GetCurrentVolume()->GetName();
}

//_____________________________________________________________________________
const char* TGeant3TGeo::CurrentVolOffName(Int_t off) const
{
  //
  // Return the current volume "off" upward in the geometrical tree
  // ID, name and copy number
  // if name=0 no name is returned
  //
  Int_t i;
  if( (i=fGcvolu->nlevel-off-1) < 0 ) {
    Warning("CurrentVolOffName",
	    "Offset requested %d but stack depth %d\n",off,fGcvolu->nlevel);
    return 0;
  }
  if (off<0 || off>gGeoManager->GetLevel()) return 0;
  if (off==0) return CurrentVolName();
  TGeoNode *node = gGeoManager->GetMother(off);
  if (!node) return 0;
  return node->GetVolume()->GetName();
}

//______________________________________________________________________
const char* TGeant3TGeo::CurrentVolPath()
{
// Return the path in geometry tree for the current volume
// ---

  return GetPath();
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::VolId(const Text_t *name) const
{
  //
  // Return the unique numeric identifier for volume name
  //
  char sname[20];
  Int_t len = strlen(name)-1;
  if (name[len] != ' ') return fMCGeo->VolId(name);
  strncpy(sname, name, len);
  sname[len] = 0;
  return fMCGeo->VolId(sname);
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::NofVolumes() const
{
  //
  // Return total number of volumes in the geometry
  //
  return fMCGeo->NofVolumes();
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::NofVolDaughters(const char* volName) const
{
// Return number of daughters of the volume specified by volName
// According to A. Morsch' G3toRoot class
// ---

  return fMCGeo->NofVolDaughters(volName);
}

//_____________________________________________________________________________
const char*  TGeant3TGeo::VolDaughterName(const char* volName, Int_t i) const
{
// Return the name of i-th daughters of the volume specified by volName
// According to A. Morsch' G3toRoot class
// ---

  return fMCGeo->VolDaughterName(volName, i);
}


//_____________________________________________________________________________
Int_t TGeant3TGeo::VolDaughterCopyNo(const char* volName, Int_t i) const
{
// Return the copyNo of i-th daughters of the volume specified by volName
// According to A. Morsch' G3toRoot class
// ---

  return fMCGeo->VolDaughterCopyNo(volName, i);
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::VolId2Mate(Int_t id) const
{
  //
  // Return material number for a given volume id
  //
  return fMCGeo->VolId2Mate(id);
}

//_____________________________________________________________________________
const char* TGeant3TGeo::VolName(Int_t id) const
{
  //
  // Return the volume name given the volume identifier
  //
  return fMCGeo->VolName(id);
}

//_____________________________________________________________________________
void TGeant3TGeo::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
		       Double_t dens, Double_t radl, Double_t absl, Float_t* buf,
		       Int_t nwbuf)
{
  //
  // Defines a Material
  //
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //

  Float_t* fbuf = CreateFloatArray(buf, nwbuf);
  G3Material(kmat, name, a, z, dens, radl, absl, fbuf, nwbuf);
  delete [] fbuf;

  fMCGeo->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
}

//_____________________________________________________________________________
void TGeant3TGeo::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
		       Double_t dens, Double_t radl, Double_t absl, Double_t* buf,
		       Int_t nwbuf)
{
  //
  // Defines a Material
  //
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //


  Float_t* fbuf = CreateFloatArray(buf, nwbuf);
  G3Material(kmat, name, a, z, dens, radl, absl, fbuf, nwbuf);
  delete [] fbuf;

  fMCGeo->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
}

//_____________________________________________________________________________
void TGeant3TGeo::Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z,
		      Double_t dens, Int_t nlmat, Float_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  //
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture.
  //
  // If nlmat < 0 then WMAT contains the number of atoms
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //

  G3Mixture(kmat, name, a, z, dens, nlmat, wmat);

  fMCGeo->Mixture(kmat, name, a, z, dens, TMath::Abs(nlmat), wmat);
}

//_____________________________________________________________________________
void TGeant3TGeo::Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z,
		      Double_t dens, Int_t nlmat, Double_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  //
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture.
  //
  // If nlmat < 0 then WMAT contains the number of atoms
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //

  Float_t* fa = CreateFloatArray(a, TMath::Abs(nlmat));
  Float_t* fz = CreateFloatArray(z, TMath::Abs(nlmat));
  Float_t* fwmat = CreateFloatArray(wmat, TMath::Abs(nlmat));

  G3Mixture(kmat, name, fa, fz, dens, nlmat, fwmat);
  Int_t i;
  for (i=0; i<TMath::Abs(nlmat); i++) {
    a[i] = fa[i]; z[i] = fz[i]; wmat[i] = fwmat[i];
  }

  delete [] fa;
  delete [] fz;
  delete [] fwmat;

  fMCGeo->Mixture(kmat, name, a, z, dens, TMath::Abs(nlmat), wmat);
}

//_____________________________________________________________________________
void TGeant3TGeo::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
		     Int_t ifield, Double_t fieldm, Double_t tmaxfd,
		     Double_t stemax, Double_t deemax, Double_t epsil,
		     Double_t stmin, Float_t* ubuf, Int_t nbuf)
{
  //
  //  kmed      tracking medium number assigned
  //  name      tracking medium name
  //  nmat      material number
  //  isvol     sensitive volume flag
  //  ifield    magnetic field
  //  fieldm    max. field value (kilogauss)
  //  tmaxfd    max. angle due to field (deg/step)
  //  stemax    max. step allowed
  //  deemax    max. fraction of energy lost in a step
  //  epsil     tracking precision (cm)
  //  stmin     min. step due to continuous processes (cm)
  //
  //  ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
  //  ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
  //  performed with g3helix; ifield = 3 if tracking performed with g3helx3.
  //

  G3Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil,
           stmin, ubuf, nbuf);

  fMCGeo->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax,
                 epsil, stmin, ubuf, nbuf);
}

//_____________________________________________________________________________
void TGeant3TGeo::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
		     Int_t ifield, Double_t fieldm, Double_t tmaxfd,
		     Double_t stemax, Double_t deemax, Double_t epsil,
		     Double_t stmin, Double_t* ubuf, Int_t nbuf)
{
  //
  //  kmed      tracking medium number assigned
  //  name      tracking medium name
  //  nmat      material number
  //  isvol     sensitive volume flag
  //  ifield    magnetic field
  //  fieldm    max. field value (kilogauss)
  //  tmaxfd    max. angle due to field (deg/step)
  //  stemax    max. step allowed
  //  deemax    max. fraction of energy lost in a step
  //  epsil     tracking precision (cm)
  //  stmin     min. step due to continuos processes (cm)
  //
  //  ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
  //  ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
  //  performed with g3helix; ifield = 3 if tracking performed with g3helx3.
  //

  Float_t* fubuf = CreateFloatArray(ubuf, nbuf);
  G3Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil,
           stmin, fubuf, nbuf);
  delete [] fubuf;

  fMCGeo->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax,
                 epsil, stmin, ubuf, nbuf);
}

//_____________________________________________________________________________
void TGeant3TGeo::Matrix(Int_t& krot, Double_t thex, Double_t phix, Double_t they,
		     Double_t phiy, Double_t thez, Double_t phiz)
{
  //
  //  krot     rotation matrix number assigned
  //  theta1   polar angle for axis i
  //  phi1     azimuthal angle for axis i
  //  theta2   polar angle for axis ii
  //  phi2     azimuthal angle for axis ii
  //  theta3   polar angle for axis iii
  //  phi3     azimuthal angle for axis iii
  //
  //  it defines the rotation matrix number irot.
  //
  krot = -1;
  fMCGeo->Matrix(krot, thex, phix, they, phiy, thez, phiz);
}

//_____________________________________________________________________________
void  TGeant3TGeo::SetRootGeometry()
{
// Notify Geant3 about use of TGeo geometry.
// The materials and tracking medias will be imported from
// TGeo at FinishGeometry().


  fImportRootGeometry = kTRUE;
}  

//_____________________________________________________________________________
const char *TGeant3TGeo::GetPath()
{
// Get current path inside G3 geometry
   return gGeoManager->GetPath();
}

//_____________________________________________________________________________
const char *TGeant3TGeo::GetNodeName()
{
// Get name of current G3 node
   if (gGeoManager->IsOutside()) return "";
   return gGeoManager->GetCurrentNode()->GetName();
}

//______________________________________________________________________
Bool_t TGeant3TGeo::GetTransformation(const TString &volumePath,TGeoHMatrix &mat)
{
    // Returns the Transformation matrix between the volume specified
    // by the path volumePath and the Top or mater volume. The format
    // of the path volumePath is as follows (assuming ALIC is the Top volume)
    // "/ALIC_1/DDIP_1/S05I_2/S05H_1/S05G_3". Here ALIC is the top most
    // or master volume which has only 1 instance of. Of all of the daughter
    // volumes of ALICE, DDIP volume copy #1 is indicated. Similarly for
    // the daughter volume of DDIP is S05I copy #2 and so on.
    // Inputs:
    //   TString& volumePath  The volume path to the specific volume
    //                        for which you want the matrix. Volume name
    //                        hierarchy is separated by "/" while the
    //                        copy number is appended using a "_".
    // Outputs:
    //  TGeoHMatrix &mat      A matrix with its values set to those
    //                        appropriate to the Local to Master transformation
    // Return:
    //   A logical value if kFALSE then an error occurred and no change to
    //   mat was made.

   // We have to preserve the modeler state
   return fMCGeo->GetTransformation(volumePath, mat);
}   
   
//______________________________________________________________________
Bool_t TGeant3TGeo::GetShape(const TString &volumePath,TString &shapeType,
                         TArrayD &par)
{
    // Returns the shape and its parameters for the volume specified
    // by volumeName.
    // Inputs:
    //   TString& volumeName  The volume name
    // Outputs:
    //   TString &shapeType   Shape type
    //   TArrayD &par         A TArrayD of parameters with all of the
    //                        parameters of the specified shape.
    // Return:
    //   A logical indicating whether there was an error in getting this
    //   information
   return fMCGeo->GetShape(volumePath, shapeType, par);
}
   
//______________________________________________________________________
Bool_t TGeant3TGeo::GetMaterial(const TString &volumeName,
                            TString &name,Int_t &imat,
                            Double_t &a,Double_t &z,Double_t &dens,
                            Double_t &radl,Double_t &inter,TArrayD &par)
{
    // Returns the Material and its parameters for the volume specified
    // by volumeName.
    // Note, Geant3 stores and uses mixtures as an element with an effective
    // Z and A. Consequently, if the parameter Z is not integer, then
    // this material represents some sort of mixture.
    // Inputs:
    //   TString& volumeName  The volume name
    // Outputs:
    //   TSrting   &name       Material name
    //   Int_t     &imat       Material index number
    //   Double_t  &a          Average Atomic mass of material
    //   Double_t  &z          Average Atomic number of material
    //   Double_t  &dens       Density of material [g/cm^3]
    //   Double_t  &radl       Average radiation length of material [cm]
    //   Double_t  &inter      Average interaction length of material [cm]
    //   TArrayD   &par        A TArrayD of user defined parameters.
    // Return:
    //   kTRUE if no errors
   Int_t i,jma,nbuf;
   Float_t af,zf,densf,radlf,interf;
   Float_t *ubuf;
   Char_t namec[20] = {20*'\0'};
   TGeoVolume *vol = gGeoManager->GetVolume(volumeName.Data());
   if (!vol) return kFALSE;
   TGeoMedium *med = vol->GetMedium();
   if (!med) return kFALSE;
   TGeoMaterial *mat = med->GetMaterial();
   imat = mat->GetUniqueID();   

   nbuf = jma = Lq()[Gclink()->jmate-imat];
   ubuf = new Float_t[nbuf];
   Gfmate(imat,namec,af,zf,densf,radlf,interf,ubuf,nbuf);
   name = mat->GetName();
   name = name.Strip(TString::kTrailing, '$');
   //
   par.Set(nbuf);
   for(i=0;i<nbuf;i++) par.AddAt(((Double_t)ubuf[i]),i);
   delete[] ubuf;
   a      = mat->GetA();
   z      = mat->GetZ();
   dens   = mat->GetDensity();
   radl   = radlf;
   inter  = interf;
   return kTRUE;
}

//______________________________________________________________________
Bool_t TGeant3TGeo::GetMedium(const TString &volumeName,TString &name,
                          Int_t &imed,Int_t &nmat,Int_t &isvol,Int_t &ifield,
                          Double_t &fieldm,Double_t &tmaxfd,Double_t &stemax,
                          Double_t &deemax,Double_t &epsil, Double_t &stmin,
                          TArrayD &par)
{
    // Returns the Medium and its parameters for the volume specified
    // by volumeName.
    // Inputs:
    //   TString& volumeName  The volume name.
    // Outputs:
    //   TString  &name       Medium name
    //   Int_t    &nmat       Material number defined for this medium
    //   Int_t    &imed       The medium index number
    //   Int_t    &isvol      volume number defined for this medium
    //   Int_t    &iflield    Magnetic field flag
    //   Double_t &fieldm     Magnetic field strength
    //   Double_t &tmaxfd     Maximum angle of deflection per step
    //   Double_t &stemax     Maximum step size
    //   Double_t &deemax     Maximum fraction of energy allowed to be lost
    //                        to continuous process.
    //   Double_t &epsil      Boundary crossing precision
    //   Double_t &stmin      Minimum step size allowed
    //   TArrayD  &par        A TArrayD of user parameters with all of the
    //                        parameters of the specified medium.
    // Return:
    //   kTRUE if there where no errors
   Int_t i,nbuf;
   Float_t fieldmf,tmaxfdf,stemaxf,deemaxf,epsilf,stminf;
   Float_t *buf;
   Char_t namec[25] = {25*'\0'};
   TGeoVolume *vol = gGeoManager->GetVolume(volumeName.Data());
   if (!vol) return kFALSE;
   TGeoMedium *med = vol->GetMedium();
   if (!med) return kFALSE;
   imed = med->GetId();
   nbuf = Lq()[Gclink()->jtmed-imed];
   buf  = new Float_t[nbuf];
   Gftmed(imed,namec,nmat,isvol,ifield,fieldmf,tmaxfdf,stemaxf,deemaxf,
          epsilf,stminf,buf,&nbuf);
   name = med->GetName();
   name = name.Strip(TString::kTrailing, '$');
   par.Set(nbuf);
   for(i=0;i<nbuf;i++) par.AddAt(((Double_t)buf[i]),i);
   delete[] buf;
   fieldm = (Double_t) fieldmf;
   tmaxfd = (Double_t) tmaxfdf;
   stemax = (Double_t) stemaxf;
   deemax = (Double_t) deemaxf;
   epsil  = (Double_t) epsilf;
   stmin  = (Double_t) stminf;
   return kTRUE;
}         

//_____________________________________________________________________________
Int_t  TGeant3TGeo::GetMedium() const
{
// Temporary added 

  return TGeant3::GetMedium();
}
  

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GBASE
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//_____________________________________________________________________________
void  TGeant3TGeo::Ggclos()
{
  //
  //   Closes off the geometry setting.
  //   Initializes the search list for the contents of each
  //   volume following the order they have been positioned, and
  //   inserting the content '0' when a call to GSNEXT (-1) has
  //   been required by the user.
  //   Performs the development of the JVOLUM structure for all
  //   volumes with variable parameters, by calling GGDVLP.
  //   Interprets the user calls to GSORD, through GGORD.
  //   Computes and stores in a bank (next to JVOLUM mother bank)
  //   the number of levels in the geometrical tree and the
  //   maximum number of contents per level, by calling GGNLEV.
  //   Sets status bit for CONCAVE volumes, through GGCAVE.
  //   Completes the JSET structure with the list of volume names
  //   which identify uniquely a given physical detector, the
  //   list of bit numbers to pack the corresponding volume copy
  //   numbers, and the generic path(s) in the JVOLUM tree,
  //   through the routine GHCLOS.
  //
  fVolNames = 0;
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gprint(const char * /*name*/)
{
  //
  // Routine to print data structures
  // CHNAME   name of a data structure
  //
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void  TGeant3TGeo::Gsmate(Int_t imat, const char *name, Float_t a, Float_t z,
		   Float_t dens, Float_t radl, Float_t absl)
{
  //
  // Defines a Material
  //
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //
  Float_t *ubuf=0;
  Int_t   nbuf=0;
  if (dens <= 0 && a != 0 && z != 0) {
     Warning("Gsmate","Density was o, set to 0.01 for imat=%d, name=%s",imat,name);
     dens = 0.01;
  }
  g3smate(imat,PASSCHARD(name), a, z, dens, radl, absl, ubuf, nbuf
	 PASSCHARL(name));

  gGeoManager->Material(name,a,z,dens,imat);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsmixt(Int_t imat, const char *name, Float_t *a, Float_t *z,
		   Float_t dens, Int_t nlmat, Float_t *wmat)
{
  //
  //       Defines mixture OR COMPOUND IMAT as composed by
  //       THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  //
  //       If NLMAT.GT.0 then WMAT contains the PROPORTION BY
  //       WEIGTHS OF EACH BASIC MATERIAL IN THE MIXTURE.
  //
  //       If NLMAT.LT.0 then WMAT contains the number of atoms
  //       of a given kind into the molecule of the COMPOUND
  //       In this case, WMAT in output is changed to relative
  //       weigths.
  //
  g3smixt(imat,PASSCHARD(name),a,z,dens,nlmat,wmat PASSCHARL(name));
  fMCGeo->Mixture(imat, name, a, z, dens, TMath::Abs(nlmat), wmat);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gstmed(Int_t numed, const char *name, Int_t nmat, Int_t isvol,
		      Int_t ifield, Float_t fieldm, Float_t tmaxfd,
		      Float_t stemax, Float_t deemax, Float_t epsil,
		      Float_t stmin)
{
  //
  //  NTMED  Tracking medium number
  //  NAME   Tracking medium name
  //  NMAT   Material number
  //  ISVOL  Sensitive volume flag
  //  IFIELD Magnetic field
  //  FIELDM Max. field value (Kilogauss)
  //  TMAXFD Max. angle due to field (deg/step)
  //  STEMAX Max. step allowed
  //  DEEMAX Max. fraction of energy lost in a step
  //  EPSIL  Tracking precision (cm)
  //  STMIN  Min. step due to continuous processes (cm)
  //
  //  IFIELD = 0 if no magnetic field; IFIELD = -1 if user decision in GUSWIM;
  //  IFIELD = 1 if tracking performed with G3RKUTA; IFIELD = 2 if tracking
  //  performed with G3HELIX; IFIELD = 3 if tracking performed with G3HELX3.
  //
  Float_t *ubuf=0;
  Int_t   nbuf=0;
  g3stmed(numed,PASSCHARD(name), nmat, isvol, ifield, fieldm, tmaxfd, stemax,
	 deemax, epsil, stmin, ubuf, nbuf PASSCHARL(name));

  gGeoManager->Medium(name,numed,nmat, isvol, ifield, fieldm, tmaxfd, stemax,deemax, epsil, stmin);
}
//_____________________________________________________________________________
void  TGeant3TGeo::Gtreve()
{
  //
  //   Controls tracking of all particles belonging to the current event
  //
  g3treve();
}

//_____________________________________________________________________________
void  TGeant3TGeo::GtreveRoot()
{
  //
  //   Controls tracking of all particles belonging to the current event
  //
  gtreveroot();
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GGEOM
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void  TGeant3TGeo::Gdtom(Float_t *xd, Float_t *xm, Int_t iflag)
{
  //
  //  Computes coordinates XM (Master Reference System
  //  knowing the coordinates XD (Detector Ref System)
  //  The local reference system can be initialized by
  //    - the tracking routines and GDTOM used in GUSTEP
  //    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
  //        (inverse routine is GMTOD)
  //
  //   If IFLAG=1  convert coordinates
  //      IFLAG=2  convert direction cosinus
  //
   Double_t XM[3], XD[3];
   Int_t i;
   for (i=0;i<3;i++) XD[i] = xd[i];
   if (iflag == 1) gGeoManager->LocalToMaster(XD,XM);
   else            gGeoManager->LocalToMasterVect(XD,XM);
   for (i=0;i<3;i++) xm[i]=XM[i];
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gdtom(Double_t *xd, Double_t *xm, Int_t iflag)
{
  //
  //  Computes coordinates XM (Master Reference System
  //  knowing the coordinates XD (Detector Ref System)
  //  The local reference system can be initialized by
  //    - the tracking routines and GDTOM used in GUSTEP
  //    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
  //        (inverse routine is GMTOD)
  //
  //   If IFLAG=1  convert coordinates
  //      IFLAG=2  convert direction cosinus
  //

   if (iflag == 1) gGeoManager->LocalToMaster(xd,xm);
   else            gGeoManager->LocalToMasterVect(xd,xm);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gmedia(Float_t *x, Int_t &numed)
{
  //
  //   Finds in which volume/medium the point X is, and updates the
  //    common /GCVOLU/ and the structure JGPAR accordingly.
  //
  //   NUMED returns the tracking medium number, or 0 if point is
  //         outside the experimental setup.
  //

   gCurrentNode = gGeoManager->FindNode(x[0],x[1],x[2]);
   if (gGeoManager->IsOutside()) {
      numed=0;
   } else {
      gcvolu->nlevel = 1 + gGeoManager->GetLevel();
      gGeoManager->GetBranchNames(gcvolu->names);
      gGeoManager->GetBranchNumbers(gcvolu->number,gcvolu->lvolum);
      TGeoVolume *vol = gCurrentNode->GetVolume();
      if (vol) {
         TGeoMedium *medium = vol->GetMedium();
         if (medium) numed = medium->GetId();
      } else {
         printf("ERROR: gmedia: NULL volume\n");
      }
   }
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gmtod(Float_t *xm, Float_t *xd, Int_t iflag)
{
  //
  //       Computes coordinates XD (in DRS)
  //       from known coordinates XM in MRS
  //       The local reference system can be initialized by
  //         - the tracking routines and GMTOD used in GUSTEP
  //         - a call to GMEDIA(XM,NUMED,CHECK)
  //         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER)
  //             (inverse routine is GDTOM)
  //
  //        If IFLAG=1  convert coordinates
  //           IFLAG=2  convert direction cosinus
  //
   Double_t XM[3], XD[3];
   Int_t i;
   for (i=0;i<3;i++) XM[i]=xm[i];
   if (iflag == 1) gGeoManager->MasterToLocal(XM,XD);
   else            gGeoManager->MasterToLocalVect(XM,XD);
   for (i=0;i<3;i++) xd[i] = XD[i];
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gmtod(Double_t *xm, Double_t *xd, Int_t iflag)
{
  //
  //       Computes coordinates XD (in DRS)
  //       from known coordinates XM in MRS
  //       The local reference system can be initialized by
  //         - the tracking routines and GMTOD used in GUSTEP
  //         - a call to GMEDIA(XM,NUMED,CHECK)
  //         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER)
  //             (inverse routine is GDTOM)
  //
  //        If IFLAG=1  convert coordinates
  //           IFLAG=2  convert direction cosinus
  //


   if (iflag == 1) gGeoManager->MasterToLocal(xm,xd);
   else            gGeoManager->MasterToLocalVect(xm,xd);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsdvn(const char *name, const char *mother, Int_t ndiv,
		     Int_t iaxis)
{
  //
  // Create a new volume by dividing an existing one
  //
  //  NAME   Volume name
  //  MOTHER Mother volume name
  //  NDIV   Number of divisions
  //  IAXIS  Axis value
  //
  //  X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
  //  It divides a previously defined volume.
  //
  fMCGeo->Gsdvn(name, mother, ndiv, iaxis);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsdvn2(const char *name, const char *mother, Int_t ndiv,
		      Int_t iaxis, Double_t c0i, Int_t numed)
{
  //
  // Create a new volume by dividing an existing one
  //
  // Divides mother into ndiv divisions called name
  // along axis iaxis starting at coordinate value c0.
  // the new volume created will be medium number numed.
  //
  fMCGeo->Gsdvn2(name, mother, ndiv, iaxis, c0i, numed);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsdvs(const char *name, const char *mother, Float_t step,
		     Int_t iaxis, Int_t numed)
{
  //
  // Create a new volume by dividing an existing one
  //
  gGeoManager->Division(name,mother,iaxis,0,0,step,numed,"s");
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsdvs2(const char *name, const char *mother, Float_t step,
		      Int_t iaxis, Float_t c0, Int_t numed)
{
  //
  // Create a new volume by dividing an existing one
  //
  gGeoManager->Division(name,mother,iaxis,0,c0,step,numed,"sx");
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsdvt(const char *name, const char *mother, Double_t step,
		     Int_t iaxis, Int_t numed, Int_t ndvmx)
{
  //
  // Create a new volume by dividing an existing one
  //
  //       Divides MOTHER into divisions called NAME along
  //       axis IAXIS in steps of STEP. If not exactly divisible
  //       will make as many as possible and will centre them
  //       with respect to the mother. Divisions will have medium
  //       number NUMED. If NUMED is 0, NUMED of MOTHER is taken.
  //       NDVMX is the expected maximum number of divisions
  //          (If 0, no protection tests are performed)
  //
  fMCGeo->Gsdvt(name, mother, step, iaxis, numed, ndvmx);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsdvt2(const char *name, const char *mother, Double_t step,
		      Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx)
{
  //
  // Create a new volume by dividing an existing one
  //
  //           Divides MOTHER into divisions called NAME along
  //            axis IAXIS starting at coordinate value C0 with step
  //            size STEP.
  //           The new volume created will have medium number NUMED.
  //           If NUMED is 0, NUMED of mother is taken.
  //           NDVMX is the expected maximum number of divisions
  //             (If 0, no protection tests are performed)
  //
  fMCGeo->Gsdvt2(name, mother, step, iaxis, c0, numed, ndvmx);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsord(const char * /*name*/, Int_t /*iax*/)
{
  //
  //    Flags volume CHNAME whose contents will have to be ordered
  //    along axis IAX, by setting the search flag to -IAX
  //           IAX = 1    X axis
  //           IAX = 2    Y axis
  //           IAX = 3    Z axis
  //           IAX = 4    Rxy (static ordering only  -> GTMEDI)
  //           IAX = 14   Rxy (also dynamic ordering -> GTNEXT)
  //           IAX = 5    Rxyz (static ordering only -> GTMEDI)
  //           IAX = 15   Rxyz (also dynamic ordering -> GTNEXT)
  //           IAX = 6    PHI   (PHI=0 => X axis)
  //           IAX = 7    THETA (THETA=0 => Z axis)
  //

}

//_____________________________________________________________________________
void  TGeant3TGeo::Gspos(const char *name, Int_t nr, const char *mother, Double_t x,
		     Double_t y, Double_t z, Int_t irot, const char *konly)
{
  //
  // Position a volume into an existing one
  //
  //  NAME   Volume name
  //  NUMBER Copy number of the volume
  //  MOTHER Mother volume name
  //  X      X coord. of the volume in mother ref. sys.
  //  Y      Y coord. of the volume in mother ref. sys.
  //  Z      Z coord. of the volume in mother ref. sys.
  //  IROT   Rotation matrix number w.r.t. mother ref. sys.
  //  ONLY   ONLY/MANY flag
  //
  //  It positions a previously defined volume in the mother.
  //
  fMCGeo->Gspos(name, nr, mother, x, y, z, irot, konly);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsposp(const char *name, Int_t nr, const char *mother,
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Float_t *upar, Int_t np )
{
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //

  fMCGeo->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsposp(const char *name, Int_t nr, const char *mother,
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Double_t *upar, Int_t np )
{
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //

 fMCGeo->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gsrotm(Int_t nmat, Float_t theta1, Float_t phi1, Float_t theta2,
		      Float_t phi2, Float_t theta3, Float_t phi3)
{
  //
  //  nmat   Rotation matrix number
  //  THETA1 Polar angle for axis I
  //  PHI1   Azimuthal angle for axis I
  //  THETA2 Polar angle for axis II
  //  PHI2   Azimuthal angle for axis II
  //  THETA3 Polar angle for axis III
  //  PHI3   Azimuthal angle for axis III
  //
  //  It defines the rotation matrix number IROT.
  //

  gGeoManager->Matrix(nmat, theta1, phi1, theta2, phi2, theta3, phi3);
}

//_____________________________________________________________________________
void  TGeant3TGeo::Gprotm(Int_t nmat)
{
  //
  //    To print rotation matrices structure JROTM
  //     nmat     Rotation matrix number
  //
  TIter next(gGeoManager->GetListOfMatrices());
  TGeoMatrix *matrix;
  while ((matrix = (TGeoMatrix*)next())) {
     if (UInt_t(nmat) == matrix->GetUniqueID()) {
        matrix->Print();
	return;
     }
  }     
  Error("Gprotm","Rotation with id=%i not found", nmat);
 }

//_____________________________________________________________________________
Int_t TGeant3TGeo::Gsvolu(const char *name, const char *shape, Int_t nmed,
		      Float_t *upar, Int_t npar)
{
  //
  //  NAME   Volume name
  //  SHAPE  Volume type
  //  NUMED  Tracking medium number
  //  NPAR   Number of shape parameters
  //  UPAR   Vector containing shape parameters
  //
  //  It creates a new volume in the JVOLUM data structure.
  //

  Int_t ivolu = 0;
  ivolu = fMCGeo->Gsvolu(name, shape, nmed, upar, npar);
  return ivolu;

}

//_____________________________________________________________________________
Int_t TGeant3TGeo::Gsvolu(const char *name, const char *shape, Int_t nmed,
		      Double_t *upar, Int_t npar)
{
  //
  //  NAME   Volume name
  //  SHAPE  Volume type
  //  NUMED  Tracking medium number
  //  NPAR   Number of shape parameters
  //  UPAR   Vector containing shape parameters
  //
  //  It creates a new volume in the JVOLUM data structure.
  //


  Int_t ivolu = 0;
  ivolu = fMCGeo->Gsvolu(name, shape, nmed, upar, npar);
  return ivolu;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//           T H E    D R A W I N G   P A C K A G E
//           ======================================
//  Drawing functions. These functions allow the visualization in several ways
//  of the volumes defined in the geometrical data structure. It is possible
//  to draw the logical tree of volumes belonging to the detector (DTREE),
//  to show their geometrical specification (DSPEC,DFSPC), to draw them
//  and their cut views (DRAW, DCUT). Moreover, it is possible to execute
//  these commands when the hidden line removal option is activated; in
//  this case, the volumes can be also either translated in the space
//  (SHIFT), or clipped by boolean operation (CVOL). In addition, it is
//  possible to fill the surfaces of the volumes
//  with solid colours when the shading option (SHAD) is activated.
//  Several tools (ZOOM, LENS) have been developed to zoom detailed parts
//  of the detectors or to scan physical events as well.
//  Finally, the command MOVE will allow the rotation, translation and zooming
//  on real time parts of the detectors or tracks and hits of a simulated event.
//  Ray-tracing commands. In case the command (DOPT RAYT ON) is executed,
//  the drawing is performed by the Geant ray-tracing;
//  automatically, the color is assigned according to the tracking medium of each
//  volume and the volumes with a density lower/equal than the air are considered
//  transparent; if the option (USER) is set (ON) (again via the command (DOPT)),
//  the user can set color and visibility for the desired volumes via the command
//  (SATT), as usual, relatively to the attributes (COLO) and (SEEN).
//  The resolution can be set via the command (SATT * FILL VALUE), where (VALUE)
//  is the ratio between the number of pixels drawn and 20 (user coordinates).
//  Parallel view and perspective view are possible (DOPT PROJ PARA/PERS); in the
//  first case, we assume that the first mother volume of the tree is a box with
//  dimensions 10000 X 10000 X 10000 cm and the view point (infinetely far) is
//  5000 cm far from the origin along the Z axis of the user coordinates; in the
//  second case, the distance between the observer and the origin of the world
//  reference system is set in cm by the command (PERSP NAME VALUE); grand-angle
//  or telescopic effects can be achieved changing the scale factors in the command
//  (DRAW). When the final picture does not occupy the full window,
//  mapping the space before tracing can speed up the drawing, but can also
//  produce less precise results; values from 1 to 4 are allowed in the command
//  (DOPT MAPP VALUE), the mapping being more precise for increasing (VALUE); for
//  (VALUE = 0) no mapping is performed (therefore max precision and lowest speed).
//  The command (VALCUT) allows the cutting of the detector by three planes
//  ortogonal to the x,y,z axis. The attribute (LSTY) can be set by the command
//  SATT for any desired volume and can assume values from 0 to 7; it determines
//  the different light processing to be performed for different materials:
//  0 = dark-matt, 1 = bright-matt, 2 = plastic, 3 = ceramic, 4 = rough-metals,
//  5 = shiny-metals, 6 = glass, 7 = mirror. The detector is assumed to be in the
//  dark, the ambient light luminosity is 0.2 for each basic hue (the saturation
//  is 0.9) and the observer is assumed to have a light source (therefore he will
//  produce parallel light in the case of parallel view and point-like-source
//  light in the case of perspective view).
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void TGeant3TGeo::Gsatt(const char *name, const char *att, Int_t val)
{
  //
  //  NAME   Volume name
  //  IOPT   Name of the attribute to be set
  //  IVAL   Value to which the attribute is to be set
  //
  //  name= "*" stands for all the volumes.
  //  iopt can be chosen among the following :
  //
  //     WORK   0=volume name is inactive for the tracking
  //            1=volume name is active for the tracking (default)
  //
  //     SEEN   0=volume name is invisible
  //            1=volume name is visible (default)
  //           -1=volume invisible with all its descendants in the tree
  //           -2=volume visible but not its descendants in the tree
  //
  //     LSTY   line style 1,2,3,... (default=1)
  //            LSTY=7 will produce a very precise approximation for
  //            revolution bodies.
  //
  //     LWID   line width -7,...,1,2,3,..7 (default=1)
  //            LWID<0 will act as abs(LWID) was set for the volume
  //            and for all the levels below it. When SHAD is 'ON', LWID
  //            represent the linewidth of the scan lines filling the surfaces
  //            (whereas the FILL value represent their number). Therefore
  //            tuning this parameter will help to obtain the desired
  //            quality/performance ratio.
  //
  //     COLO   colour code -166,...,1,2,..166 (default=1)
  //            n=1=black
  //            n=2=red;    n=17+m, m=0,25, increasing luminosity according to 'm';
  //            n=3=green;  n=67+m, m=0,25, increasing luminosity according to 'm';
  //            n=4=blue;   n=117+m, m=0,25, increasing luminosity according to 'm';
  //            n=5=yellow; n=42+m, m=0,25, increasing luminosity according to 'm';
  //            n=6=violet; n=142+m, m=0,25, increasing luminosity according to 'm';
  //            n=7=lightblue; n=92+m, m=0,25, increasing luminosity according to 'm';
  //            colour=n*10+m, m=1,2,...9, will produce the same colour
  //            as 'n', but with increasing luminosity according to 'm';
  //            COLO<0 will act as if abs(COLO) was set for the volume
  //            and for all the levels below it.
  //            When for a volume the attribute FILL is > 1 (and the
  //            option SHAD is on), the ABS of its colour code must be < 8
  //            because an automatic shading of its faces will be
  //            performed.
  //
  //     FILL  (1992) fill area  -7,...,0,1,...7 (default=0)
  //            when option SHAD is "on" the FILL attribute of any
  //            volume can be set different from 0 (normal drawing);
  //            if it is set to 1, the faces of such volume will be filled
  //            with solid colours; if ABS(FILL) is > 1, then a light
  //            source is placed along the observer line, and the faces of
  //            such volumes will be painted by colours whose luminosity
  //            will depend on the amount of light reflected;
  //            if ABS(FILL) = 1, then it is possible to use all the 166
  //            colours of the colour table, becouse the automatic shading
  //            is not performed;
  //            for increasing values of FILL the drawing will be performed
  //            with higher and higher resolution improving the quality (the
  //            number of scan lines used to fill the faces increases with FILL);
  //            it is possible to set different values of FILL
  //            for different volumes, in order to optimize at the same time
  //            the performance and the quality of the picture;
  //            FILL<0 will act as if abs(FILL) was set for the volume
  //            and for all the levels below it.
  //            This kind of drawing can be saved in 'picture files'
  //            or in view banks.
  //            0=drawing without fill area
  //            1=faces filled with solid colours and resolution = 6
  //            2=lowest resolution (very fast)
  //            3=default resolution
  //            4=.................
  //            5=.................
  //            6=.................
  //            7=max resolution
  //            Finally, if a coloured background is desired, the FILL
  //            attribute for the first volume of the tree must be set
  //            equal to -abs(colo), colo being >0 and <166.
  //
  //     SET   set number associated to volume name
  //     DET   detector number associated to volume name
  //     DTYP  detector type (1,2)
  //
//  InitHIGZ();
  gGeoManager->SetVolumeAttribute(name, att, val);
}

//_____________________________________________________________________________
Int_t TGeant3TGeo::Glvolu(Int_t /*nlev*/, Int_t * /*lnam*/,Int_t * /*lnum*/)
{
  //
  //  nlev   number of leveles deap into the volume tree
  //         size of the arrays lnam and lnum
  //  lnam   an integer array whos 4 bytes contain the askii code for the
  //         volume names
  //  lnum   an integer array containing the copy numbers for that specific
  //         volume
  //
  //  This routine fills the volulme paramters in common /gcvolu/ for a
  //  physical tree, specified by the list lnam and lnum of volume names
  //  and numbers, and for all its ascendants up to level 1. This routine
  //  is optimsed and does not re-compute the part of the history already
  //  available in GCVOLU. This means that if it is used in user programs
  //  outside the usual framwork of the tracking, the user has to initilise
  //  to zero NLEVEL in the common GCVOLU. It return 0 if there were no
  //  problems in make the call.
  //
  return 0;
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdshow(Int_t /*iview*/)
{
  //
  //  IVIEW  View number
  //
  //  It shows on the screen the contents of a view bank. It
  //  can be called after a view bank has been closed.
  //
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdopt(const char * /*name*/,const char * /*value*/)
{
  //
  //  NAME   Option name
  //  VALUE  Option value
  //
  //  To set/modify the drawing options.
  //     IOPT   IVAL      Action
  //
  //     THRZ    ON       Draw tracks in R vs Z
  //             OFF (D)  Draw tracks in X,Y,Z
  //             180
  //             360
  //     PROJ    PARA (D) Parallel projection
  //             PERS     Perspective
  //     TRAK    LINE (D) Trajectory drawn with lines
  //             POIN       " " with markers
  //     HIDE    ON       Hidden line removal using the CG package
  //             OFF (D)  No hidden line removal
  //     SHAD    ON       Fill area and shading of surfaces.
  //             OFF (D)  Normal hidden line removal.
  //     RAYT    ON       Ray-tracing on.
  //             OFF (D)  Ray-tracing off.
  //     EDGE    OFF      Does not draw contours when shad is on.
  //             ON  (D)  Normal shading.
  //     MAPP    1,2,3,4  Mapping before ray-tracing.
  //             0   (D)  No mapping.
  //     USER    ON       User graphics options in the raytracing.
  //             OFF (D)  Automatic graphics options.
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdraw(const char * /*name*/,Double_t /*theta*/, Double_t /*phi*/, Double_t /*psi*/,
			Double_t /*u0*/,Double_t /*v0*/,Double_t /*ul*/,Double_t /*vl*/)
{
  //
  //  NAME   Volume name
  //  +
  //  THETA  Viewing angle theta (for 3D projection)
  //  PHI    Viewing angle phi (for 3D projection)
  //  PSI    Viewing angle psi (for 2D rotation)
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  This function will draw the volumes,
  //  selected with their graphical attributes, set by the Gsatt
  //  facility. The drawing may be performed with hidden line removal
  //  and with shading effects according to the value of the options HIDE
  //  and SHAD; if the option SHAD is ON, the contour's edges can be
  //  drawn or not. If the option HIDE is ON, the detector can be
  //  exploded (BOMB), clipped with different shapes (CVOL), and some
  //  of its parts can be shifted from their original
  //  position (SHIFT). When HIDE is ON, if
  //  the drawing requires more than the available memory, the program
  //  will evaluate and display the number of missing words
  //  (so that the user can increase the
  //  size of its ZEBRA store). Finally, at the end of each drawing (with HIDE on),
  //  the program will print messages about the memory used and
  //  statistics on the volumes' visibility.
  //  The following commands will produce the drawing of a green
  //  volume, specified by NAME, without using the hidden line removal
  //  technique, using the hidden line removal technique,
  //  with different linewidth and colour (red), with
  //  solid colour, with shading of surfaces, and without edges.
  //  Finally, some examples are given for the ray-tracing. (A possible
  //  string for the NAME of the volume can be found using the command DTREE).
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdrawc(const char * /*name*/,Int_t /*axis*/, Float_t /*cut*/,Float_t /*u0*/,
			 Float_t /*v0*/,Float_t /*ul*/,Float_t /*vl*/)
{
  //
  //  NAME   Volume name
  //  CAXIS  Axis value
  //  CUTVAL Cut plane distance from the origin along the axis
  //  +
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  The cut plane is normal to caxis (X,Y,Z), corresponding to iaxis (1,2,3),
  //  and placed at the distance cutval from the origin.
  //  The resulting picture is seen from the the same axis.
  //  When HIDE Mode is ON, it is possible to get the same effect with
  //  the CVOL/BOX function.
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdrawx(const char * /*name*/,Float_t /*cutthe*/, Float_t /*cutphi*/,
		     Float_t /*cutval*/, Float_t /*theta*/, Float_t /*phi*/, Float_t /*u0*/,
			 Float_t /*v0*/,Float_t /*ul*/,Float_t /*vl*/)
{
  //
  //  NAME   Volume name
  //  CUTTHE Theta angle of the line normal to cut plane
  //  CUTPHI Phi angle of the line normal to cut plane
  //  CUTVAL Cut plane distance from the origin along the axis
  //  +
  //  THETA  Viewing angle theta (for 3D projection)
  //  PHI    Viewing angle phi (for 3D projection)
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  The cut plane is normal to the line given by the cut angles
  //  cutthe and cutphi and placed at the distance cutval from the origin.
  //  The resulting picture is seen from the viewing angles theta,phi.
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdspec(const char * /*name*/)
{
  //
  //  NAME   Volume name
  //
  //  Shows 3 views of the volume (two cut-views and a 3D view), together with
  //  its geometrical specifications. The 3D drawing will
  //  be performed according the current values of the options HIDE and
  //  SHAD and according the current SetClipBox clipping parameters for that
  //  volume.
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::DrawOneSpec(const char * /*name*/)
{
  //
  //  Function called when one double-clicks on a volume name
  //  in a TPavelabel drawn by Gdtree.
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::Gdtree(const char * /*name*/,Int_t /*levmax*/, Int_t /*isel*/)
{
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree,
  //  Each volume in the tree is represented by a TPaveTree object.
  //  Double-clicking on a TPaveTree draws the specs of the corresponding volume.
  //  Use TPaveTree pop-up menu to select:
  //    - drawing specs
  //    - drawing tree
  //    - drawing tree of parent
  //
  return;
}

//_____________________________________________________________________________
void TGeant3TGeo::GdtreeParent(const char * /*name*/,Int_t /*levmax*/, Int_t /*isel*/)
{
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree of the parent of name.
  //
}

//____________________________________________________________________________
Int_t  TGeant3TGeo::ImportMaterial(const TGeoMaterial* mat)
{
// Imports the Root material in Geant3 and returns its Geant3 index
// ---

  Int_t kmat;
  const TGeoMixture* mixt = dynamic_cast<const TGeoMixture*>(mat);
  if (mixt) {
    // TGeo stores only proportions by weigth
    Int_t nlmat = mixt->GetNelements();
    Float_t* fa = CreateFloatArray(mixt->GetAmixt(), TMath::Abs(nlmat));
    Float_t* fz = CreateFloatArray(mixt->GetZmixt(), TMath::Abs(nlmat));
    Float_t* fwmat = CreateFloatArray(mixt->GetWmixt(), TMath::Abs(nlmat));
    G3Mixture(kmat, mixt->GetName(), fa, fz, mixt->GetDensity(), TMath::Abs(nlmat), fwmat);
    delete [] fa;
    delete [] fz;
    delete [] fwmat;
  }
  else {
    Float_t* buf = 0;
    // Inject radlen with negative sign to be stored in G3
    Double_t radlen = mat->GetRadLen();
    // Ignore abslen from TGeo and let G3 compute it
    G3Material(kmat, mat->GetName(), mat->GetA(), mat->GetZ(),
               mat->GetDensity(), -radlen, 0, buf, 0);
  }
  return kmat;
}

//____________________________________________________________________________
void TGeant3TGeo::FinishGeometry()
{
  //
  // Finalise geometry construction
  //

  //Close the geometry structure
  if (gDebug > 0) printf("FinishGeometry, calling ggclos\n");
  Ggclos();

  if (fImportRootGeometry) {

    // Import materials
    //
    TIter next1(gGeoManager->GetListOfMaterials());
    TGeoMaterial* mat;
    Int_t nofMaterials = 0;
    while ((mat=(TGeoMaterial*)next1())) {
      Int_t kmat = ImportMaterial(mat);
      mat->SetUniqueID(kmat);
      nofMaterials++;
    }  	         

    // Number of media
    Int_t nofMedia = 0;
    TIter next2(gGeoManager->GetListOfMedia());
    TGeoMedium* medx;
    while ((medx=(TGeoMedium*)next2())) nofMedia++;

    // Import media
    //
    Int_t  maxNofMaterials = nofMaterials + nofMedia;
    TArrayI usedMaterials(maxNofMaterials);
    for (Int_t i=0; i<maxNofMaterials; i++)
      usedMaterials[i] = -1;

    TIter next3(gGeoManager->GetListOfMedia());
    TGeoMedium* med;
    while ((med=(TGeoMedium*)next3())) {
      Int_t kmed;
      Int_t nmat = med->GetMaterial()->GetUniqueID();

      // if material is already used define a new Geant3 material
      // (do not reset TGeoMaterial index)
      if (usedMaterials[nmat] >0 )
        nmat = ImportMaterial(med->GetMaterial());
      usedMaterials[nmat] = 1;        

      Int_t isvol  = (Int_t) med->GetParam(0);
      Int_t ifield = (Int_t) med->GetParam(1);
      Double_t fieldm = med->GetParam(2);
      Double_t tmaxfd = med->GetParam(3);
      Double_t stemax = med->GetParam(4);
      Double_t deemax = med->GetParam(5);
      Double_t epsil  = med->GetParam(6);
      Double_t stmin  = med->GetParam(7);
      G3Medium(kmed, med->GetName(), nmat, isvol, ifield, fieldm, tmaxfd,
               stemax,deemax, epsil, stmin);
      med->SetId(kmed);
    }
    if (gDebug > 0) printf("FinishGeometry, geometry retreived from file, materials/media mapped to G3\n");
  } else {
    TGeoVolume *top = (TGeoVolume*)gGeoManager->GetListOfVolumes()->First();
    gGeoManager->SetTopVolume(top);
    if (gDebug > 0) printf("FinishGeometry, calling CloseGeometry\n");
    gGeoManager->CloseGeometry();
  }

  if (gDebug > 0) printf("FinishGeometry, calling MisalignGeometry()\n");
  TVirtualMCApplication::Instance()->MisalignGeometry();
  //  gROOT->GetListOfBrowsables()->Add(gGeoManager);
  if (gDebug > 0) printf("FinishGeometry, calling SetColors\n");

  //Create the color table
  SetColors();
  if (gDebug > 0) printf("FinishGeometry, returning\n");
}

//_____________________________________________________________________________
void TGeant3TGeo::SetColors()
{
  //
  // Set the colors for all the volumes
  // this is done sequentially for all volumes
  // based on the number of their medium
  //
  TIter next(gGeoManager->GetListOfVolumes());
  TGeoVolume *volume;
  while ((volume = (TGeoVolume*)next())) {
     if (volume->IsAssembly()) continue;
     TGeoMedium *medium = (TGeoMedium*)volume->GetMedium();
     Int_t icol = medium->GetId()%6+2;
     volume->SetLineColor(icol);
  }
}

//_____________________________________________________________________________
//
//                 Interfaces to Fortran
//
//_____________________________________________________________________________


//______________________________________________________________________
void ginvolTGeo(Float_t *x, Int_t &isame)
{
   if (gGeoManager->IsSameLocation(x[0], x[1], x[2])) isame = 1;
   else isame = 0;
}


//______________________________________________________________________
void gtmediTGeo(Float_t *x, Int_t &numed)
{
   gcchan->lsamvl = kTRUE;
   gCurrentNode = gGeoManager->FindNode(x[0],x[1],x[2]);
   gcchan->lsamvl = gGeoManager->IsSameLocation();
   if (gGeoManager->IsOutside()) {
      numed=0;
   } else {
      gcvolu->nlevel = 1 + gGeoManager->GetLevel();
      gGeoManager->GetBranchNames(gcvolu->names);
      gGeoManager->GetBranchNumbers(gcvolu->number,gcvolu->lvolum);
      TGeoVolume *vol = gCurrentNode->GetVolume();
      if (vol) {
         TGeoMedium *medium = vol->GetMedium();
         if (medium) numed = medium->GetId();
      } else {
         printf("ERROR: gtmedi: NULL volume\n");
      }
   }
}


//______________________________________________________________________
void gmediaTGeo(Float_t *x, Int_t &numed, Int_t & /*check*/)
{
   gCurrentNode = gGeoManager->FindNode(x[0],x[1],x[2]);
   if (gGeoManager->IsOutside()) {
      numed=0;
   } else {
      gcvolu->nlevel = 1 + gGeoManager->GetLevel();
      gGeoManager->GetBranchNames(gcvolu->names);
      gGeoManager->GetBranchNumbers(gcvolu->number,gcvolu->lvolum);
      TGeoVolume *vol = gCurrentNode->GetVolume();
      if (vol) {
         TGeoMedium *medium = vol->GetMedium();
         if (medium) numed = medium->GetId();
      } else {
         printf("ERROR: gmedia: NULL volume\n");
      }
   }
}

//______________________________________________________________________
void gtmanyTGeo(Int_t &level1)
{
   if (level1==1) {
      Int_t nlevel = gcvolu->nlevel;
//      printf("gtmanyTGeo nlevel=%i %s\n",nlevel,gGeoManager->GetPath());
      gcvol1->nlevl1 = nlevel;
      if (nlevel>0) {
         memcpy(gcvol1->names1, gcvolu->names, nlevel*sizeof(Int_t));
         memcpy(gcvol1->numbr1, gcvolu->number, nlevel*sizeof(Int_t));
         memcpy(gcvol1->lvolu1, gcvolu->lvolum, nlevel*sizeof(Int_t));
      }  
   }   
}

//______________________________________________________________________
void gtonlyTGeo(Int_t &isOnly)
{
   //with Geant3, return gonly(nlevel);
//   if (gGeoManager->IsCurrentOverlapping()) isOnly = 0;
//   else isOnly = 1;
   // With TGeo, G3 is seeing a ONLY geometry
   isOnly = 1;
}

//_____________________________________________________________________________
void glvoluTGeo(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier)
{
  //
  //  nlev   number of levels deap into the volume tree
  //         size of the arrays lnam and lnum
  //  lnam   an integer array whos 4 bytes contain the askii code for the
  //         volume names
  //  lnum   an integer array containing the copy numbers for that specific
  //         volume
  //
  //  This routine fills the volume parameters in common /gcvolu/ for a
  //  physical tree, specified by the list lnam and lnum of volume names
  //  and numbers, and for all its ascendants up to level 1. This routine
  //  is optimised and does not re-compute the part of the history already
  //  available in GCVOLU. This means that if it is used in user programs
  //  outside the usual framework of the tracking, the user has to initialize
  //  to zero NLEVEL in the common GCVOLU. It return 0 if there were no
  //  problems in make the call.
  //
   TGeoVolume *vol = gGeoManager->GetTopVolume();
   TGeoVolume *vdaughter = 0;
   TGeoNode *node = 0;
   Int_t nd;
   Bool_t found = kFALSE;
   ier = 0;
   gGeoManager->CdTop();
   if (nlev<1) nlev = 1;
   gcvolu->nlevel = nlev;
   if (nlev==1) return;
   Int_t *lvol = gcvol1->lvolu1;
   memcpy(gcvolu->names, lnam,  nlev*sizeof(Int_t));
   memcpy(gcvolu->number, lnum, nlev*sizeof(Int_t));
   memcpy(gcvolu->lvolum, lvol, nlev*sizeof(Int_t));
//   for (Int_t i=0;i<nlev;i++) printf(" #%i: %i  %i  %i\n", i, lnam[i], lnum[i],lvol[i]);
   
   for (Int_t i=1; i<nlev; i++) {
      nd = vol->GetNdaughters();
      found = kFALSE;
      for (Int_t id=0; id<nd; id++) {
         node = vol->GetNode(id);
         vdaughter = node->GetVolume();
         if (vdaughter->GetNumber() == lvol[i]) {
            if (node->GetNumber()==lnum[i]) {
               found = kTRUE;
               gGeoManager->CdDown(id);
               vol = vdaughter;
               break;
            } 
         }
      }
      if (!found) {
         printf("### ERROR in TGeant3TGeo::glvoluTGeo(): cannot restore path\n");
         ier = 1;
         return;
      }           
   }   
}


//______________________________________________________________________
void gtnextTGeo()
{
   Float_t *x = gctrak->vect;
   Double_t step = gctrak->step;
   Int_t itrtyp = gckine->itrtyp;
   gGeoManager->SetCurrentPoint(x[0],x[1],x[2]);
   gGeoManager->SetCurrentDirection(x[3],x[4],x[5]);
   if (step<=0) {
      gctrak->safety = 0.;
      gctrak->snext = 0.;
      gctrak->ignext = 0;
      return;
   }
   // Find distance to next boundary. Global matrix computed only if
   // gtnext is called by gtckov.
   if (itrtyp==7) gGeoManager->FindNextBoundary(-step);
   else           gGeoManager->FindNextBoundary(step);
   gctrak->safety = gGeoManager->GetSafeDistance();
   Double_t snext  = gGeoManager->GetStep();
   if (snext<=0) {
      gctrak->safety = 0.;
      gctrak->snext = 0.;
      gctrak->ignext = 1;
      return;
   }
   if (snext < step) {
      gctrak->snext  = snext;
      gctrak->ignext = 1;
   } else {
      gctrak->ignext = 0;
      gctrak->snext = gctrak->step;
   }
}


//______________________________________________________________________
void ggperpTGeo(Float_t * /*x*/, Float_t *norm, Int_t &ierr)
{
// Computes the normal to the next crossed surface, assuming that
// FindNextBoundary() was already called.
   ierr = 0;
   Double_t *dblnorm = gGeoManager->FindNormalFast();
   if (!dblnorm) {
      ierr = 1;
      return;
   }
   norm[0] = -dblnorm[0];
   norm[1] = -dblnorm[1];
   norm[2] = -dblnorm[2];
}

