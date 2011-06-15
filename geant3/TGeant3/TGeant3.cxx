
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

/* $Id$ */

/*
$Log: TGeant3.cxx,v $
Revision 1.59  2007/07/25 20:06:49  brun
From Ivana:
Changed the fPDGCode type to TArrayI to allow
its dynamic size (thanks to Susan Kasahara).

Revision 1.58  2007/07/24 19:43:24  brun
From Ivana:
Do not add particles in TDatabasePDG if they
are were already aded;
(thanks to Susan Kasahara for this suggestion)

Revision 1.57  2007/07/23 20:04:03  brun
From Ivana:
Implemented the new pdg "standard" codes for ions
defined in pdg-2006, thanks to Susan Kasahara
for this suggestion.

Revision 1.56  2007/05/18 08:44:15  brun
A major update of GEANTE by Andrea Fontana and Alberto Rotondi


1) update of the Coulomb multiple scattering parametrization;
2) update of the straggling of energy loss for thin materials;
3) new options to extrapolate the track parameters to the point
   of closest approach to a point or to a wire (straight line).

Details on the physical motivation behind this work can be found
in our report for the Panda Collaboration, available at:

http://www.pv.infn.it/~fontana/tracking.pdf

Feel free to contact us for questions and discussions about these
topics by using the following email addresses:

alberto.rotondi@pv.infn.it
andrea.fontana@pv.infn.it

---

List of changes in the fortran and C++ routines of the geant3
VMC directory:

- gcmore.inc
  gtmore.inc
  geant3LinkDef.h
  gcomad.F

 Added a new common that contains all the new variables:
      COMMON/GCMORE/GCALPHA,ICLOSE,PFINAL(3),DSTRT,WIRE1(3),WIRE2(3),
     +              P1(3),P2(3),P3(3),CLENG(3)

     input to ERLAND:
      GCALPHA: energy cut parameter for energy loss fluctuations

     input to EUSTEP:
      ICLOSE: = 1 the use of the common is enabled for the closest
                  approach to a point PFINAL(3)
              = 2 the use of the common is enabled for the closest
                  approach to a wire of extremes WIRE1(3) and WIRE2(3)
              = 0 the common is empty and disabled
      PFINAL(3): assigned point
      DSTRT: assigned distance between initial point in ERTRAK
             and PFINAL along straight line (currently noy used)
      WIRE1(3): first point of a wire
      WIRE2(3): second point of a wire

     output from EUSTEP:
      P1(3): point previous to the point of closest approach to
             PFINAL() or wire
      P2(3): point of closest approach to PFINAL() or wire
      P3(3): point next to the point of closest approach to
             PFINAL() or wire
      CLENG(3): track length to the previous 3 points

      Important note: the calculated points of closest approach are
      depending on the GEANE steps. For calculating the true point
      of closest approach the last 3 points of the extrapolation, i.e.
      the previous to closest, the closest and the next to closest are
      returned to the user. Different algorithms can be implemented, but
      we decided to leave this to the users in the C++ interface to GEANE.

- ermcsc.F
 new expression for the variance of the Coulomb multiple scattering
 according to Fruhwirth and Regler, NIM A 456 (2001) 369

- ertrch.F
 added DESTEP in the calling string of ERLAND for calculation with
 Urban model. Added and saved previous step PRSTEP.

- erland.F
 added new calculation for sigma of straggling in energy loss
 to include in Geane the Urban/Landau approximation, as explained
 in the Geant manual and related papers.
 The model parametrization can be controlled with a user variable (GCALPHA)
 in the new GCMORE common block: 1 is for old gaussian model valid
 for dense materials, other values (see the report) are for gaseous
 materials.

- eustep.F
 added the calculation to the distance of closest approach to a point
 or to a wire.

- TGeant3.h
- TGeant3.cxx
 added the possibility to define user cuts (already present in the gccuts
 struct but not in the TGeant3::SetCUTS method) and to define the new
 variables of the GCMORE common with two new methods SetECut() and
 SetClose().
 Added new method InitGEANE() to initialize GEANE to the old behaviour
 (default) for backward compatibility. Only the multiple scattering has
 been updated to a more correct formula.
 Corrected a typo in the call to the routine Trscsd().

Revision 1.55  2007/03/26 10:15:04  brun
Fix a problem when adding a new tracking medium to the TObjArray.

Revision 1.54  2007/03/23 21:11:44  brun
From Ivana Hrivnacova and Andrea Fontana:
Reintroduce functionality in TGeant3::SetCuts that was removed in a previous patch.

Revision 1.53  2007/03/22 08:58:41  brun
From Ivana:
Restore the function TGeant3::MediumId

Revision 1.52  2007/02/28 16:25:14  brun
From Federico:
Suppress compiler warnings coming from unused arguments

Revision 1.51  2006/12/19 13:16:19  brun
from Mohammad Al-Turany & Denis Bertini

Changes in  TGeant3/TGeant3.cxx and TGeant3.h
------------------------------------
1. Geane interface functions are added:
    void  eufill(Int_t n, Float_t *ein, Float_t *xlf);
    void  eufilp(const int n,Float_t *ein,Float_t *pli,Float_t *plf);
    void  eufilv(Int_t n, Float_t *ein,	Char_t *namv, Int_t *numv,Int_t *iovl);
    void  trscsd(Float_t *pc,Float_t *rc,Float_t *pd,Float_t *rd,
		 Float_t *h,Float_t ch,Int_t ierr,Float_t spu,Float_t *dj,Float_t *dk);
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

    Ertrio_t  *fErtrio
    Ertrio1_t *fErtrio1
    Eropts_t  *fEropts
    Eroptc_t  *fEroptc
    Erwork_t  *fErwork
    Trcom3_t  *fTrcom3

2. The size of the error matrix errin is corrected to 15

Revision 1.49  2006/04/20 10:14:55  brun
From Peter Hristov:
small change in TGeant3.cxx:
to avoid some rare problems in glandz.F and probably in other places where the random number should not
be equal to 0 and 1.

Revision 1.48  2006/03/15 08:12:11  brun
-New Makefile.macosx

-delete all TGeant3 GUI and graphics classes

-Remove all references to GUI/Graphics classes from G3Volume.cxx and TGeant3.cxx

-Remove TGeant3/galicef.F and rdummies.F

-Remove references to GUI/Graphics classes from geant3LinkDef.h

-delete added/dummies2.F

-change added/dummies.c, keeping only __attribute__ and MAIN__

-Change Makefile accordingly

Revision 1.47  2005/11/18 21:25:22  brun
From Bjorn, Andrei:
Implemented new VMC functions for access to geometry;
added -Woverloaded-virtual to Makefile.linux

Revision 1.46  2005/07/21 17:54:37  brun
Implement same code in the float* versions that were already implemented
in the Double* versions.

Revision 1.45  2005/07/20 09:22:50  brun
From Federico:
Fixes to compile with gcc4CVS: ----------------------------------------------------------------------

Revision 1.44  2005/05/21 05:48:33  brun
In TGeant3::Init put the code calling ConstructGeometry conditional
on ROOT versions >-5.01/01

Revision 1.43  2005/05/17 12:47:00  brun
From Ivana:
- Added call to new TVirtualMCApplication::ConstructOpGeometry() function
- Bug fix in G3Mixture (do not update wmat values twice)
- In CreateFloatArray: added check for values > FLT_MAX

Revision 1.42  2005/02/16 08:01:51  brun
Implement function grndmq. It returns gRandom->GetSeed().
The result makes sense only when Trandom is the generator.

Revision 1.41  2005/02/08 11:22:03  brun
From Ivana:
For TGeant3.h:
Added IsRootGeometrySupported() function
(now required by TVirtualMC)

For TGeant3.cxx:
Updated text in Fatal in SetRootGeometry.

Revision 1.40  2005/01/07 11:34:43  brun
Change the J/Psi lifetime from 0 to 7.6e-21
(thanks Yuri Kharlov)

Revision 1.39  2004/12/21 15:34:48  brun
Implement TGeant3TGeo::isRootGeometry returning kTRUE

Revision 1.38  2004/12/17 11:55:47  brun
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
Add a protection in TGeant3::Gsmate in case the material density is null

Revision 1.16  2003/12/01 23:51:22  brun
>From Andrei and Peter:
add a few missing cases when compiling with the WITHROOT option.

Revision 1.15  2003/11/28 09:44:15  brun
New version of TGeant3 supporting the options WITHG3 and WITHROOT

Revision 1.14  2003/10/09 06:28:45  brun
In TGeant3::ParticleName, increase size of local array name[20] to name[21]

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
    TMCParticleType ParticleMCType(..) const
- corrected charge in AddParticlesToPdgDataBase

Revision 1.12  2003/07/22 06:53:28  brun
This version does not yet support TGeo geometry.
TVirtualMC must be initialized with the 3rd argument set to kFALSE

Revision 1.11  2003/07/18 10:22:50  brun
Changes to reflect the equivalent changes in the abstract classes in vmc
(thanks Peter Hristov)

Revision 1.10  2003/07/16 07:40:09  brun
>From Andreas Morsch

- default g3 specific initialisation moved to TGeant3::Init()
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

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Interface Class to the Geant3.21 Monte Carlo                      //
//                                                                    //
//                                                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <ctype.h>
#include <stdlib.h>

#include "TROOT.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TString.h"
#include "TParameter.h"
#include "TGeoMatrix.h"
#include "TObjString.h"

#include "TGeant3.h"

#include "TCallf77.h"
#include "TVirtualMCDecayer.h"
#include "TPDGCode.h"

#ifndef WIN32
# define g3zebra  g3zebra_
# define grfile   grfile_
# define g3pcxyz  g3pcxyz_
# define g3gclos  g3gclos_
# define g3last   g3last_
# define g3init   g3init_
# define g3cinit  g3cinit_
# define g3run    g3run_
# define g3trig   g3trig_
# define g3trigc  g3trigc_
# define g3trigi  g3trigi_
# define g3work   g3work_
# define g3zinit  g3zinit_
# define g3fmate  g3fmate_
# define g3fpart  g3fpart_
# define g3ftmed  g3ftmed_
# define g3ftmat  g3ftmat_
# define g3mate   g3mate_
# define g3part   g3part_
# define g3sdk    g3sdk_
# define g3smate  g3smate_
# define g3fang   g3fang_ 
# define g3smixt  g3smixt_
# define g3spart  g3spart_
# define g3stmed  g3stmed_
# define g3sckov  g3sckov_
# define g3stpar  g3stpar_
# define g3fkine  g3fkine_
# define g3fvert  g3fvert_
# define g3skine  g3skine_
# define g3svert  g3svert_
# define g3physi  g3physi_
# define g3debug  g3debug_
# define g3ekbin  g3ekbin_
# define g3finds  g3finds_
# define g3sking  g3sking_
# define g3skpho  g3skpho_
# define g3sstak  g3sstak_
# define g3sxyz   g3sxyz_
# define g3many   g3many_
# define g3track  g3track_
# define g3treve  g3treve_
# define gtreveroot  gtreveroot_
# define grndm    grndm_
# define grndmq   grndmq_
# define g3dtom   g3dtom_
# define g3lmoth  g3lmoth_
# define g3media  g3media_
# define g3mtod   g3mtod_
# define g3sdvn   g3sdvn_
# define g3sdvn2  g3sdvn2_
# define g3sdvs   g3sdvs_
# define g3sdvs2  g3sdvs2_
# define g3sdvt   g3sdvt_
# define g3sdvt2  g3sdvt2_
# define g3sord   g3sord_
# define g3spos   g3spos_
# define g3sposp  g3sposp_
# define g3srotm  g3srotm_
# define g3protm  g3protm_
# define g3svolu  g3svolu_
# define g3print  g3print_
# define dzshow   dzshow_
# define g3satt   g3satt_
# define g3fpara  g3fpara_
# define gckpar   gckpar_
# define g3ckmat  g3ckmat_
# define g3lvolu  g3lvolu_
# define geditv   geditv_
# define mzdrop   mzdrop_

# define ertrak  ertrak_
# define ertrgo  ertrgo_
# define eufill  eufill_
# define eufilp  eufilp_
# define eufilv  eufilv_
# define trscsp  trscsp_
# define trspsc  trspsc_
# define trscsd  trscsd_
# define trsdsc  trsdsc_
# define erxyzc  erxyzc_

# define gcomad gcomad_

# define g3brelm g3brelm_
# define g3prelm g3prelm_

# define rxgtrak rxgtrak_
# define rxouth  rxouth_
# define rxinh   rxinh_


#else

# define gzebra  GZEBRA
# define grfile  GRFILE
# define gpcxyz  GPCXYZ
# define ggclos  GGCLOS
# define glast   GLAST
# define ginit   GINIT
# define g3cinit G3CINIT
# define grun    GRUN
# define gtrig   GTRIG
# define gtrigc  GTRIGC
# define gtrigi  GTRIGI
# define gwork   GWORK
# define g3zinit G3ZINIT
# define gfmate  GFMATE
# define gfpart  GFPART
# define gftmed  GFTMED
# define gftmat  GFTMAT
# define gmate   GMATE
# define gpart   GPART
# define gsdk    GSDK
# define gsmate  GSMATE
# define gsmixt  GSMIXT
# define gspart  GSPART
# define gstmed  GSTMED
# define gsckov  GSCKOV
# define gstpar  GSTPAR
# define gfkine  GFKINE
# define gfvert  GFVERT
# define gskine  GSKINE
# define gsvert  GSVERT
# define gphysi  GPHYSI
# define gdebug  GDEBUG
# define gekbin  GEKBIN
# define gfinds  GFINDS
# define gsking  GSKING
# define gskpho  GSKPHO
# define gsstak  GSSTAK
# define gsxyz   GSXYZ
# define gtrack  GTRACK
# define gtreve  GTREVE
# define gtreveroot  GTREVEROOT
# define grndm   GRNDM
# define grndmq  GRNDMQ
# define gdtom   GDTOM
# define glmoth  GLMOTH
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
# define dzshow  DZSHOW
# define gsatt   GSATT
# define gfpara  GFPARA
# define gckpar  GCKPAR
# define gckmat  GCKMAT
# define glvolu  GLVOLU
# define geditv  GEDITV
# define mzdrop  MZDROP

# define ertrak  ERTRAK
# define ertrgo  ERTRGO
# define eufill  EUFILL
# define eufilp  EUFILP
# define eufilv  EUFILV
# define trscsp  TRSCSP
# define trspsc  TRSPSC
# define trscsd  TRSCSD
# define trsdsc  TRSDSC
# define erxyzc  ERXYZC

# define gcomad  GCOMAD

# define gbrelm  GBRELM
# define gprelm  GPRELM

# define rxgtrak RXGTRAK
# define rxouth  RXOUTH
# define rxinh   RXINH
# define gfang   GFANG 

#endif

//______________________________________________________________________
extern "C"
{
  //
  // Prototypes for GEANT functions
  //
  void type_of_call g3zebra(const int&);

  void type_of_call g3pcxyz();

  void type_of_call g3gclos();

  void type_of_call g3last();

  void type_of_call g3init();

  void type_of_call g3cinit();

  void type_of_call g3run();

  void type_of_call g3trig();

  void type_of_call g3trigc();

  void type_of_call g3trigi();

  void type_of_call g3work(const int&);

  void type_of_call g3zinit();

  void type_of_call g3mate();

  void type_of_call g3part();

  void type_of_call g3sdk(Int_t &, Float_t *, Int_t *);

  void type_of_call g3fkine(Int_t &, Float_t *, Float_t *, Int_t &,
			   Int_t &, Float_t *, Int_t &);

  void type_of_call g3fvert(Int_t &, Float_t *, Int_t &, Int_t &,
			   Float_t &, Float_t *, Int_t &);

  void type_of_call g3skine(Float_t *,Int_t &, Int_t &, Float_t *,
			   Int_t &, Int_t &);

  void type_of_call g3svert(Float_t *,Int_t &, Int_t &, Float_t *,
			   Int_t &, Int_t &);

  void type_of_call g3physi();

  void type_of_call g3debug();

  void type_of_call g3ekbin();

  void type_of_call g3finds();

  void type_of_call g3sking(Int_t &);

  void type_of_call g3skpho(Int_t &);

  void type_of_call g3sstak(Int_t &);

  void type_of_call g3sxyz();

  void type_of_call g3many();

  void type_of_call g3track();

  void type_of_call g3treve();

  void type_of_call gtreveroot();

  void type_of_call grndm(Float_t *r, const Int_t &n)
  {
     //gRandom->RndmArray(n,r);
     for(Int_t i=0; i<n; i++)
       do r[i]=gRandom->Rndm(); while(0>=r[i] || r[i]>=1);
  }

  void type_of_call grndmq(Int_t &is1, Int_t &is2, const Int_t &,
			   DEFCHARD DEFCHARL)
  {is1=gRandom->GetSeed(); is2=0; /*only valid with TRandom;*/}

  void type_of_call g3dtom(Float_t *, Float_t *, Int_t &);

  void type_of_call g3lmoth(DEFCHARD, Int_t &, Int_t &, Int_t *,
			   Int_t *, Int_t * DEFCHARL);

  void type_of_call g3media(Float_t *, Int_t &, Int_t&);

  void type_of_call g3mtod(Float_t *, Float_t *, Int_t &);

  void type_of_call g3srotm(const Int_t &, const Float_t &, const Float_t &,
			   const Float_t &, const Float_t &, const Float_t &,
			   const Float_t &);

  void type_of_call g3protm(const Int_t &);

  void type_of_call g3rfile(const Int_t&, DEFCHARD,
			   DEFCHARD DEFCHARL DEFCHARL);

  void type_of_call g3fmate(const Int_t&, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t &, Float_t *,
			   Int_t& DEFCHARL);

  void type_of_call g3fang( Float_t *, Float_t &,
			    Float_t &, Float_t &, Float_t &,
			    Int_t & ); 

  void type_of_call g3fpart(const Int_t&, DEFCHARD, Int_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL);

  void type_of_call g3ftmed(const Int_t&, DEFCHARD, Int_t &, Int_t &, Int_t &,
			   Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t * DEFCHARL);

  void type_of_call g3ftmat(const Int_t&, const Int_t&, DEFCHARD, const Int_t&,
			   Float_t*, Float_t*
			   ,Float_t *, Int_t & DEFCHARL);

  void type_of_call g3smate(const Int_t&, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t &, Float_t *,
			   Int_t & DEFCHARL);

  void type_of_call g3smixt(const Int_t&, DEFCHARD, const Float_t *,
               const Float_t *, const Float_t &, const Int_t &,
               Float_t * DEFCHARL);

  void type_of_call g3spart(const Int_t&, DEFCHARD, Int_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL);


  void type_of_call g3stmed(const Int_t&, DEFCHARD, Int_t &, Int_t &, Int_t &,
			   Float_t &, Float_t &, Float_t &, Float_t &,
			   Float_t &, Float_t &, Float_t *, Int_t & DEFCHARL);

  void type_of_call g3sckov(Int_t &itmed, Int_t &npckov, Float_t *ppckov,
			   Float_t *absco, Float_t *effic, Float_t *rindex);
  void type_of_call g3stpar(const Int_t&, DEFCHARD, Float_t & DEFCHARL);

  void type_of_call g3sdvn(DEFCHARD,DEFCHARD, Int_t &, Int_t &
			  DEFCHARL DEFCHARL);

  void type_of_call g3sdvn2(DEFCHARD,DEFCHARD, Int_t &, Int_t &, Float_t &,
			   Int_t & DEFCHARL DEFCHARL);

  void type_of_call g3sdvs(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Int_t &
			  DEFCHARL DEFCHARL);

  void type_of_call g3sdvs2(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Float_t &,
			   Int_t & DEFCHARL DEFCHARL);

  void type_of_call g3sdvt(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Int_t &,
			  Int_t & DEFCHARL DEFCHARL);

  void type_of_call g3sdvt2(DEFCHARD,DEFCHARD, Float_t &, Int_t &, Float_t&,
			   Int_t &, Int_t & DEFCHARL DEFCHARL);

  void type_of_call g3sord(DEFCHARD, Int_t & DEFCHARL);

  void type_of_call g3spos(DEFCHARD, Int_t &, DEFCHARD, Float_t &, Float_t &,
			  Float_t &, Int_t &, DEFCHARD DEFCHARL DEFCHARL
			  DEFCHARL);

  void type_of_call g3sposp(DEFCHARD, Int_t &, DEFCHARD, Float_t &, Float_t &,
			   Float_t &, Int_t &, DEFCHARD,
			   Float_t *, Int_t & DEFCHARL DEFCHARL DEFCHARL);

  void type_of_call g3svolu(DEFCHARD, DEFCHARD, Int_t &, Float_t *, Int_t &,
			   Int_t & DEFCHARL DEFCHARL);

  void type_of_call g3satt(DEFCHARD, DEFCHARD, Int_t & DEFCHARL DEFCHARL);

  void type_of_call g3fpara(DEFCHARD , Int_t&, Int_t&, Int_t&, Int_t&, 
                            Float_t*, Float_t* DEFCHARL);

  void type_of_call gckpar(Int_t&, Int_t&, Float_t*);

  void type_of_call g3ckmat(Int_t&, DEFCHARD DEFCHARL);

  void type_of_call g3lvolu(Int_t&, Int_t*, Int_t*, Int_t&);

  void type_of_call g3print(DEFCHARD,const int& DEFCHARL);

  void type_of_call dzshow(DEFCHARD,const int&,const int&,DEFCHARD,const int&,
			   const int&, const int&, const int& DEFCHARL
			   DEFCHARL);
  void type_of_call mzdrop(Int_t&, Int_t&, DEFCHARD DEFCHARL);

  void type_of_call setbomb(Float_t &);

  void type_of_call setclip(DEFCHARD, Float_t &,Float_t &,Float_t &,Float_t &,
			    Float_t &, Float_t & DEFCHARL);

  void type_of_call gcomad(DEFCHARD, Int_t*& DEFCHARL);

  void type_of_call ertrak(const Float_t *const x1, const Float_t *const p1,
			   const Float_t *x2, const Float_t *p2,
			   const Int_t &ipa, DEFCHARD DEFCHARL);
  void type_of_call eufill(Int_t n, Float_t *ein,
		       Float_t *xlf);
  void type_of_call eufilp(const int n,Float_t *ein,
			   Float_t *pli,Float_t *plf);
  void type_of_call eufilv(Int_t n, Float_t *ein,
			Char_t *namv, Int_t *numv,Int_t *iovl);
  void type_of_call trscsd(Float_t *pc,Float_t *rc,Float_t *pd,Float_t *rd,
                         Float_t *h,Float_t ch,Int_t ierr,Float_t spu,Float_t *dj,Float_t *dk);
  void type_of_call trsdsc(Float_t *pd,Float_t *rd,Float_t *pc,Float_t *rc,
                         Float_t *h,Float_t *ch,Int_t *ierr,Float_t *spu,Float_t *dj,Float_t *dk);
  void type_of_call trscsp(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,
			Float_t *ch,Int_t *ierr,Float_t *spx);
  void type_of_call trspsc(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,
			Float_t *ch,Int_t *ierr,Float_t *spx);

  void type_of_call ertrgo();
 
  void type_of_call erxyzc();

    float type_of_call g3brelm(const Float_t &z, const Float_t& t, 
                               const Float_t& cut);
    float type_of_call g3prelm(const Float_t &z, const Float_t& t, 
                               const Float_t& cut);
}

#ifndef WIN32
#  define gudigi gudigi_
#  define guhadr guhadr_
#  define guout  guout_
#  define guphad guphad_
#  define gudcay gudcay_
#  define guiget guiget_
#  define guinme guinme_
#  define guinti guinti_
#  define gunear gunear_
#  define guskip guskip_
#  define guview guview_
#  define gupara gupara_
#  define gudtim gudtim_
#  define guplsh guplsh_
#  define gutrev gutrev_
#  define gutrak gutrak_
#  define guswim guswim_
#  define gufld  gufld_
#  define gustep gustep_
#  define gukine gukine_

#  define gheish gheish_
#  define flufin flufin_
#  define gfmfin gfmfin_
#  define gpghei gpghei_
#  define fldist fldist_
#  define gfmdis gfmdis_
#  define g3helx3 g3helx3_
#  define g3helix g3helix_
#  define g3rkuta g3rkuta_
#  define g3track g3track_
#  define gtreveroot gtreveroot_
#  define g3last  g3last_
#  define g3invol g3invol_
#  define g3tmedi g3tmedi_
#  define g3media g3media_
#  define g3tmany g3tmany_
#  define g3tnext g3tnext_
#  define g3gperp g3gperp_
#  define ginvol ginvol_
#  define gtmedi gtmedi_
#  define gtmany gtmany_
#  define gtonly gtonly_
#  define gmedia gmedia_
#  define glvolu glvolu_
#  define gtnext gtnext_
#  define ggperp ggperp_

#else
#  define gudigi GUDIGI
#  define guhadr GUHADR
#  define guout  GUOUT
#  define guphad GUPHAD
#  define gudcay GUDCAY
#  define guiget GUIGET
#  define guinme GUINME
#  define guinti GUINTI
#  define gunear GUNEAR
#  define guskip GUSKIP
#  define guview GUVIEW
#  define gupara GUPARA
#  define gudtim GUDTIM
#  define guplsh GUPLSH
#  define gutrev GUTREV
#  define gutrak GUTRAK
#  define guswim GUSWIM
#  define gufld  GUFLD
#  define gustep GUSTEP
#  define gukine GUKINE

#  define eustep EUSTEP

#  define gheish GHEISH
#  define flufin FLUFIN
#  define gfmfin GFMFIN
#  define gpghei GPGHEI
#  define fldist FLDIST
#  define gfmdis GFMDIS
#  define g3helx3 G3HELX3
#  define g3helix G3HELIX
#  define g3gperp G3GPERP
#  define g3rkuta G3RKUTA
#  define gtrack GTRACK
#  define gtreveroot GTREVEROOT
#  define glast  GLAST
#  define ginvol GINVOL
#  define gtmedi GTMEDI
#  define gtmany GTMANY
#  define gmedia GMEDIA
#  define glvolu GLVOLU
#  define gtnext GTNEXT
#  define ggperp GGPERP

#endif

extern "C" type_of_call void gheish();
extern "C" type_of_call void flufin();
extern "C" type_of_call void gfmfin();
extern "C" type_of_call void gpghei();
extern "C" type_of_call void fldist();
extern "C" type_of_call void gfmdis();
extern "C" type_of_call void g3helx3(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void g3helix(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void g3rkuta(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void g3gperp(Float_t*, Float_t*, Int_t&);
extern "C" type_of_call void g3track();
extern "C" type_of_call void gtreveroot();
extern "C" type_of_call void g3last();
extern "C" type_of_call void g3invol(Float_t*, Int_t&);
extern "C" type_of_call void g3tmedi(Float_t*, Int_t&);
extern "C" type_of_call void g3tmany(Int_t&);
extern "C" type_of_call void g3media(Float_t*, Int_t&, Int_t&);
extern "C" type_of_call void g3tnext();
extern "C" type_of_call void ginvol(Float_t*, Int_t&);
extern "C" type_of_call void gtmedi(Float_t*, Int_t&);
extern "C" type_of_call void gtmany(Int_t&);
extern "C" type_of_call void gtonly(Int_t&);
extern "C" type_of_call void gmedia(Float_t*, Int_t&, Int_t&);
extern "C" type_of_call void glvolu(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier);
extern "C" type_of_call void gtnext();
extern "C" type_of_call void ggperp(Float_t*, Float_t*, Int_t&);


//
// Geant3 global pointer
//
Gctrak_t *gctrak = 0;
Gcvolu_t *gcvolu = 0;
Gckine_t *gckine = 0;
TGeant3* geant3 = 0;
static const Int_t kDefSize = 600;
Int_t count_ginvol = 0;
Int_t count_gmedia = 0;
Int_t count_gtmedi = 0;
Int_t count_gtnext = 0;
Gcchan_t *gcchan = 0;

Gconst_t *gconst=0;          //! GCONST common structure
Gconsx_t *cconsx=0;          //! GCONSX common structure
Gcjump_t *gcjump=0;          //! GCJUMP common structure




extern "C"  type_of_call void gtonlyg3(Int_t&);
void (*fginvol)(Float_t*, Int_t&) = 0;
void (*fgtmedi)(Float_t*, Int_t&) = 0;
void (*fgtmany)(Int_t&) = 0;
void (*fgtonly)(Int_t&) = 0;
void (*fgmedia)(Float_t*, Int_t&, Int_t&) = 0;
void (*fglvolu)(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier) = 0;
void (*fgtnext)() = 0;
void (*fggperp)(Float_t*, Float_t*, Int_t&) = 0;

//#define STATISTICS
#ifdef STATISTICS
#include "TTree.h"
#include "TFile.h"
Double_t oldvect[6], oldstep, oldsafety;
Int_t statcode, statsame;
Char_t statpath[120];
Double_t statsafety, statsnext;
TTree *stattree =0;
TFile *statfile=0;
#endif

//______________________________________________________________________
TGeant3::TGeant3()
  : TVirtualMC(),
    fNG3Particles(0),
    fNPDGCodes(0),
    fPDGCode(),
    fMCGeo(0),
    fImportRootGeometry(kFALSE),
    fStopRun(kFALSE)
{
  //
  // Default constructor
  //
   geant3 = this;
}

//______________________________________________________________________
TGeant3::TGeant3(const char *title, Int_t nwgeant)
  : TVirtualMC("TGeant3",title, kFALSE),
    fNG3Particles(0),
    fNPDGCodes(0),
    fPDGCode(),
    fMCGeo(0),
    fImportRootGeometry(kFALSE),
    fStopRun(kFALSE),
    fSkipNeutrinos(kTRUE)
{
  //
  // Standard constructor for TGeant3 with ZEBRA initialization
  //

#ifdef STATISTICS
   statfile = new TFile("stat.root","recreate");
   stattree = new TTree("stat","stat tree");
   stattree->Branch("statcode",&statcode,"statcode/I");
   stattree->Branch("statsame",&statsame,"statsame/I");
   stattree->Branch("statpath",statpath,"statpath/C");
   stattree->Branch("oldvect",oldvect,"oldvect[6]/D");
   stattree->Branch("oldsafety",&oldsafety,"oldsafety/D");
   stattree->Branch("oldstep",&oldstep,"oldstep/D");
   stattree->Branch("snext",&statsnext,"statsnext/D");
   stattree->Branch("safety",&statsafety,"statsafety/D");
#endif

  geant3 = this;

  if(nwgeant) {
    g3zebra(nwgeant);
    g3init();
    g3zinit();
  } else {
    g3cinit();
  }
  //
  // Load Address of Geant3 commons
  LoadAddress();
  //
  // Zero number of particles
  fNG3Particles = 0;
  fNPDGCodes=0;
  
  // Set initial size to fPDGCode table
  fPDGCode.Set(100);

  //set pointers to tracker functions
  fginvol = g3invol;
  fgtmedi = g3tmedi;
  fgtmany = g3tmany;
  fgtonly = gtonlyg3;
  fgmedia = g3media;
  fglvolu = g3lvolu;
  fgtnext = g3tnext;
  fggperp = g3gperp;

  InitGEANE();
}

//______________________________________________________________________
TGeant3::~TGeant3()
{
  if(fVolNames) {
    delete [] fVolNames;
    fVolNames=0;
  }
}

//______________________________________________________________________
Int_t TGeant3::CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens,
			       Float_t &radl, Float_t &absl) const
{
  //
  // Return the parameters of the current material during transport
  //
  z     = fGcmate->z;
  a     = fGcmate->a;
  dens  = fGcmate->dens;
  radl  = fGcmate->radl;
  absl  = fGcmate->absl;
  return 1;  //this could be the number of elements in mixture
}

//______________________________________________________________________
void TGeant3::DefaultRange()
{
  //
  // Set range of current drawing pad to 20x20 cm
  //
}

//______________________________________________________________________
void TGeant3::InitHIGZ()
{
  //
  // Initialize HIGZ
  //
}

//______________________________________________________________________
void TGeant3::InitGEANE()
{
  //
  // Initialize GEANE for default use
  //
  Float_t pf[3]={0.,0.,0.};
  Float_t w1[3]={0.,0.,0.};
  Float_t w2[3]={0.,0.,0.};
  Float_t p1[3]={0.,0.,0.};
  Float_t p2[3]={0.,0.,0.};
  Float_t p3[3]={0.,0.,0.};
  Float_t cl[3]={0.,0.,0.};
  geant3 = this;
  geant3->SetECut(1.);
  geant3->SetClose(0,pf,999.,w1,w2,p1,p2,p3,cl);
}

//______________________________________________________________________
void TGeant3::LoadAddress()
{
  //
  // Assigns the address of the GEANT common blocks to the structures
  // that allow their access from C++
  //
   Int_t *addr;
   gcomad(PASSCHARD("QUEST"), (int*&) fQuest PASSCHARL("QUEST"));
   gcomad(PASSCHARD("GCBANK"),(int*&) fGcbank  PASSCHARL("GCBANK"));
   gcomad(PASSCHARD("GCLINK"),(int*&) fGclink  PASSCHARL("GCLINK"));
   gcomad(PASSCHARD("GCCUTS"),(int*&) fGccuts  PASSCHARL("GCCUTS"));
   gcomad(PASSCHARD("GCMORE"),(int*&) fGcmore  PASSCHARL("GCMORE"));
   gcomad(PASSCHARD("GCMULO"),(int*&) fGcmulo  PASSCHARL("GCMULO"));
   gcomad(PASSCHARD("GCFLAG"),(int*&) fGcflag  PASSCHARL("GCFLAG"));
   gcomad(PASSCHARD("GCKINE"),(int*&) fGckine  PASSCHARL("GCKINE"));
   gcomad(PASSCHARD("GCKING"),(int*&) fGcking  PASSCHARL("GCKING"));
   gcomad(PASSCHARD("GCKIN2"),(int*&) fGckin2  PASSCHARL("GCKIN2"));
   gcomad(PASSCHARD("GCKIN3"),(int*&) fGckin3  PASSCHARL("GCKIN3"));
   gcomad(PASSCHARD("GCMATE"),(int*&) fGcmate  PASSCHARL("GCMATE"));
   gcomad(PASSCHARD("GCTMED"),(int*&) fGctmed  PASSCHARL("GCTMED"));
   gcomad(PASSCHARD("GCTRAK"),(int*&) fGctrak  PASSCHARL("GCTRAK"));
   gcomad(PASSCHARD("GCTPOL"),(int*&) fGctpol  PASSCHARL("GCTPOL"));
   gcomad(PASSCHARD("GCVOLU"),(int*&) fGcvolu  PASSCHARL("GCVOLU"));
   gcomad(PASSCHARD("GCNUM"), (int*&) fGcnum   PASSCHARL("GCNUM"));
   gcomad(PASSCHARD("GCSETS"),(int*&) fGcsets  PASSCHARL("GCSETS"));
   gcomad(PASSCHARD("GCPHYS"),(int*&) fGcphys  PASSCHARL("GCPHYS"));
   gcomad(PASSCHARD("GCPHLT"),(int*&) fGcphlt  PASSCHARL("GCPHLT"));
   gcomad(PASSCHARD("GCOPTI"),(int*&) fGcopti  PASSCHARL("GCOPTI"));
   gcomad(PASSCHARD("GCTLIT"),(int*&) fGctlit  PASSCHARL("GCTLIT"));
   gcomad(PASSCHARD("GCVDMA"),(int*&) fGcvdma  PASSCHARL("GCVDMA"));
   gcomad(PASSCHARD("GCCHAN"),(int*&) gcchan   PASSCHARL("GCCHAN"));

   // Commons for GEANE
   gcomad(PASSCHARD("ERTRIO"), (int*&) fErtrio  PASSCHARL("ERTRIO"));
   gcomad(PASSCHARD("ERTRIO1"),(int*&) fErtrio1 PASSCHARL("ERTRIO1"));
   gcomad(PASSCHARD("EROPTS"), (int*&) fEropts  PASSCHARL("EROPTS"));
   gcomad(PASSCHARD("EROPTC"), (int*&) fEroptc  PASSCHARL("EROPTC"));
   gcomad(PASSCHARD("ERWORK"), (int*&) fErwork  PASSCHARL("ERWORK"));
   gcomad(PASSCHARD("GCONST"),(int*&) fGconst  PASSCHARL("GCONST"));
   gcomad(PASSCHARD("GCONSX"),(int*&) fGconsx  PASSCHARL("GCONSX"));
   gcomad(PASSCHARD("GCJUMP"),(int*&) fGcjump  PASSCHARL("GCJUMP"));

   // Variables for ZEBRA store
   gcomad(PASSCHARD("IQ"), addr  PASSCHARL("IQ"));
   fZiq = addr;
   gcomad(PASSCHARD("LQ"), addr  PASSCHARL("LQ"));
   fZlq = addr;
   fZq       = (float*)fZiq;
   gctrak = fGctrak;
   gcvolu = fGcvolu;
   gckine = fGckine;
}

//______________________________________________________________________
void TGeant3::GeomIter()
{
  //
  // Geometry iterator for moving upward in the geometry tree
  // Initialize the iterator
  //
  fNextVol=fGcvolu->nlevel;
}

//______________________________________________________________________
Int_t TGeant3::NextVolUp(Text_t *name, Int_t &copy)
{
  //
  // Geometry iterator for moving upward in the geometry tree
  // Return next volume up
  //
  fNextVol--;
  Int_t i, gname;
  if(fNextVol>=0) {
    gname=fGcvolu->names[fNextVol];
    copy=fGcvolu->number[fNextVol];
    i=fGcvolu->lvolum[fNextVol];
    name = fVolNames[i-1];
    if(gname == fZiq[fGclink->jvolum+i]) return i;
    else printf("GeomTree: Volume %s not found in bank\n",name);
  }
  return 0;
}

//______________________________________________________________________
void TGeant3::BuildPhysics()
{
  Gphysi();
}

//______________________________________________________________________
void TGeant3::AddParticlesToPdgDataBase() const
{

//
// Add particles to the PDG data base

    TDatabasePDG *pdgDB = TDatabasePDG::Instance();

    const Double_t kAu2Gev=0.9314943228;
    const Double_t khSlash = 1.0545726663e-27;
    const Double_t kErg2Gev = 1/1.6021773349e-3;
    const Double_t khShGev = khSlash*kErg2Gev;
    const Double_t kYear2Sec = 3600*24*365.25;
//
// Bottom mesons
// mass and life-time from PDG
//
// Done by default now from Pythia6 table!
//
//
// Ions
//

  if ( !pdgDB->GetParticle(GetIonPdg(1,2)) )
    pdgDB->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,
                       0,3,"Ion",GetIonPdg(1,2));

  if ( !pdgDB->GetParticle(GetIonPdg(1,3)) )
    pdgDB->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE,
                       khShGev/(12.33*kYear2Sec),3,"Ion",GetIonPdg(1,3));

  if ( !pdgDB->GetParticle(GetIonPdg(2,4)) )
    pdgDB->AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE,
                       khShGev/(12.33*kYear2Sec),6,"Ion",GetIonPdg(2,4));

  if ( !pdgDB->GetParticle(GetIonPdg(2,3)) )
    pdgDB->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,
                       0,6,"Ion",GetIonPdg(2,3));

// Special particles
//
  if ( !pdgDB->GetParticle(GetSpecialPdg(50)) )
    pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,
                       0,0,"Special",GetSpecialPdg(50));

  if ( !pdgDB->GetParticle(GetSpecialPdg(51)) )
    pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton",0,kFALSE,
                       0,0,"Special",GetSpecialPdg(51));

}


//______________________________________________________________________
Int_t TGeant3::CurrentVolID(Int_t &copy) const
{
  //
  // Returns the current volume ID and copy number
  //
  Int_t i, gname;
  if( (i=fGcvolu->nlevel-1) < 0 ) {
    Warning("CurrentVolID","Stack depth only %d\n",fGcvolu->nlevel);
  } else {
    gname=fGcvolu->names[i];
    copy=fGcvolu->number[i];
    i=fGcvolu->lvolum[i];
    if(gname == fZiq[fGclink->jvolum+i]) return i;
    else Warning("CurrentVolID","Volume %4s not found\n",(char*)&gname);
  }
  return 0;
}

//______________________________________________________________________
Int_t TGeant3::CurrentVolOffID(Int_t off, Int_t &copy) const
{
  //
  // Return the current volume "off" upward in the geometrical tree
  // ID and copy number
  //
  Int_t i, gname;
  if( (i=fGcvolu->nlevel-off-1) < 0 ) {
    Warning("CurrentVolOffID","Offset requested %d but stack depth %d\n",
	    off,fGcvolu->nlevel);
  } else {
    gname=fGcvolu->names[i];
    copy=fGcvolu->number[i];
    i=fGcvolu->lvolum[i];
    if(gname == fZiq[fGclink->jvolum+i]) return i;
    else Warning("CurrentVolOffID","Volume %4s not found\n",(char*)&gname);
  }
  return 0;
}

//______________________________________________________________________
const char* TGeant3::CurrentVolName() const
{
  //
  // Returns the current volume name
  //
  Int_t i;
  if( (i=fGcvolu->nlevel-1) < 0 ) {
    Warning("CurrentVolName","Stack depth %d\n",fGcvolu->nlevel);
    return 0;
  }
  Int_t gname=fGcvolu->names[i];
  i=fGcvolu->lvolum[i];
  if(gname == fZiq[fGclink->jvolum+i]) return fVolNames[i-1];
  else Warning("CurrentVolName","Volume %4s not found\n",(char*) &gname);
  return 0;
}

//______________________________________________________________________
const char* TGeant3::CurrentVolOffName(Int_t off) const
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
  Int_t gname=fGcvolu->names[i];
  i=fGcvolu->lvolum[i];
  if(gname == fZiq[fGclink->jvolum+i]) return fVolNames[i-1];
  else Warning("CurrentVolOffName","Volume %4s not found\n",(char*)&gname);
  return 0;
}

//______________________________________________________________________
const char* TGeant3::CurrentVolPath()
{
// Return the path in geometry tree for the current volume
// ---

  return GetPath();
}

//______________________________________________________________________
Int_t TGeant3::IdFromPDG(Int_t pdg) const
{
  //
  // Return Geant3 code from PDG and pseudo ENDF code
  //
  for(Int_t i=0;i<fNPDGCodes;++i)
    if(pdg==fPDGCode[i]) return i;
  return -1;
}

//______________________________________________________________________
Int_t TGeant3::PDGFromId(Int_t id) const
{
  //
  // Return PDG code and pseudo ENDF code from Geant3 code
  //
  if(id>0 && id<fNPDGCodes) return fPDGCode[id];
  else return -1;
}

//______________________________________________________________________
void TGeant3::DefineParticles()
{
  //
  // Define standard Geant 3 particles
  Gpart();
  //
  // Load standard numbers for GEANT particles and PDG conversion
  fPDGCode[fNPDGCodes++]=-99;   //  0 = unused location
  fPDGCode[fNPDGCodes++]=22;    //  1 = photon
  fPDGCode[fNPDGCodes++]=-11;   //  2 = positron
  fPDGCode[fNPDGCodes++]=11;    //  3 = electron
  fPDGCode[fNPDGCodes++]=12;    //  4 = neutrino e
  fPDGCode[fNPDGCodes++]=-13;   //  5 = muon +
  fPDGCode[fNPDGCodes++]=13;    //  6 = muon -
  fPDGCode[fNPDGCodes++]=111;   //  7 = pi0
  fPDGCode[fNPDGCodes++]=211;   //  8 = pi+
  fPDGCode[fNPDGCodes++]=-211;  //  9 = pi-
  fPDGCode[fNPDGCodes++]=130;   // 10 = Kaon Long
  fPDGCode[fNPDGCodes++]=321;   // 11 = Kaon +
  fPDGCode[fNPDGCodes++]=-321;  // 12 = Kaon -
  fPDGCode[fNPDGCodes++]=2112;  // 13 = Neutron
  fPDGCode[fNPDGCodes++]=2212;  // 14 = Proton
  fPDGCode[fNPDGCodes++]=-2212; // 15 = Anti Proton
  fPDGCode[fNPDGCodes++]=310;   // 16 = Kaon Short
  fPDGCode[fNPDGCodes++]=221;   // 17 = Eta
  fPDGCode[fNPDGCodes++]=3122;  // 18 = Lambda
  fPDGCode[fNPDGCodes++]=3222;  // 19 = Sigma +
  fPDGCode[fNPDGCodes++]=3212;  // 20 = Sigma 0
  fPDGCode[fNPDGCodes++]=3112;  // 21 = Sigma -
  fPDGCode[fNPDGCodes++]=3322;  // 22 = Xi0
  fPDGCode[fNPDGCodes++]=3312;  // 23 = Xi-
  fPDGCode[fNPDGCodes++]=3334;  // 24 = Omega-
  fPDGCode[fNPDGCodes++]=-2112; // 25 = Anti Neutron
  fPDGCode[fNPDGCodes++]=-3122; // 26 = Anti Lambda
  fPDGCode[fNPDGCodes++]=-3222; // 27 = Anti Sigma -
  fPDGCode[fNPDGCodes++]=-3212; // 28 = Anti Sigma 0
  fPDGCode[fNPDGCodes++]=-3112; // 29 = Anti Sigma +
  fPDGCode[fNPDGCodes++]=-3322; // 30 = Anti Xi 0      
  fPDGCode[fNPDGCodes++]=-3312; // 31 = Anti Xi +
  fPDGCode[fNPDGCodes++]=-3334; // 32 = Anti Omega +


  Int_t mode[6];
  Int_t kz, ipa;
  Float_t bratio[6];

  fNG3Particles = 33;

  /* --- Define additional particles */
  Gspart(fNG3Particles++,"OMEGA(782)",3,0.782,0.,7.836e-23);// 33 = OMEGA(782)
  fPDGCode[fNPDGCodes++]=223;   // 33 = Omega(782)

  Gspart(fNG3Particles++,"PHI(1020)",3,1.019,0.,1.486e-22);// 34 = PHI(1020)
  fPDGCode[fNPDGCodes++]=333;   // 34 = PHI (1020)

  Gspart(fNG3Particles++, "D +", 4, 1.8693, 1., 1.040e-12);      // 35 = D+        // G4DMesonPlus
  fPDGCode[fNPDGCodes++]=411;   // 35 = D+

  Gspart(fNG3Particles++, "D -", 4, 1.8693, -1., 1.040e-12);     // 36 = D-        // G4DMesonMinus
  fPDGCode[fNPDGCodes++]=-411;  // 36 = D-

  Gspart(fNG3Particles++, "D 0", 3, 1.8645, 0., 0.415e-12);      // 37 = D0        // G4DMesonZero
  fPDGCode[fNPDGCodes++]=421;   // 37 = D0

  Gspart(fNG3Particles++,"ANTI D 0",3,1.8645,0.,0.415e-12);      // 38 = Anti D0   // G4AntiDMesonZero
  fPDGCode[fNPDGCodes++]=-421;  // 38 = D0 bar


  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=-99;  // 39 = unassigned

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=-99;  // 40 = unassigned

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=-99;  // 41 = unassigned

  Gspart(fNG3Particles++, "RHO +", 4, 0.768, 1., 4.353e-24);  // 42 = Rho+
  fPDGCode[fNPDGCodes++]=213;   // 42 = RHO+

  Gspart(fNG3Particles++, "RHO -", 4, 0.768, -1., 4.353e-24); // 43 = Rho-
  fPDGCode[fNPDGCodes++]=-213;   // 43 = RHO-

  Gspart(fNG3Particles++, "RHO 0", 3, 0.768, 0., 4.353e-24);  // 44 = Rho0
  fPDGCode[fNPDGCodes++]=113;   //  44 = RHO0

//
// Ions

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=GetIonPdg(1, 2); // 45 = Deuteron

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=GetIonPdg(1, 3); // 46 = Triton

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=GetIonPdg(2, 4); // 47 = Alpha

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=0;               // 48 = geantino mapped to rootino

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=GetIonPdg(2, 3); // 49 = HE3

  fNG3Particles++;
  fPDGCode[fNPDGCodes++]=GetSpecialPdg(50); // 50 = Cherenkov
// special
  Gspart(fNG3Particles++, "FeedbackPhoton", 7, 0., 0.,1.e20 );
  fPDGCode[fNPDGCodes++]=GetSpecialPdg(51); // 51 = FeedbackPhoton
//

  Gspart(fNG3Particles++, "Lambda_c-", 4, 2.28646, +1., 2.06e-13);
  // Gspart(fNG3Particles++, "Lambda_c+", 4, 2.28646, +1., 0.200e-12); // G4LambdacPlus
  fPDGCode[fNPDGCodes++]=4122;         //52 = Lambda_c+

  Gspart(fNG3Particles++, "Lambda_c-", 4, 2.28646, -1., 2.06e-13);
  // Gspart(fNG3Particles++, "Lambda_c-", 4, 2.2849, -1., 0.200e-12);  // G4AntiLamdacPlus
  fPDGCode[fNPDGCodes++]=-4122;        //53 = Lambda_c-

  Gspart(fNG3Particles++, "D_s+", 4, 1.9682, +1., 0.490e-12);       // G4DsMesonPlus * booklet (July 2006): lifetime=0.500e-12  
  fPDGCode[fNPDGCodes++]=431;          //54 = D_s+

  Gspart(fNG3Particles++, "D_s-", 4, 1.9682, -1., 0.490e-12);       // G4DsMesonMinus * booklet: lifetime=0.500e-12
  fPDGCode[fNPDGCodes++]=-431;         //55 = D_s-

  Gspart(fNG3Particles++, "Tau+", 5, 1.77699, +1., 290.6e-15);      // G4TauPlus *
  fPDGCode[fNPDGCodes++]=-15;          //56 = Tau+

  Gspart(fNG3Particles++, "Tau-", 5, 1.77699, -1., 290.6e-15);      // G4TauMinus *
  fPDGCode[fNPDGCodes++]= 15;          //57 = Tau-

  Gspart(fNG3Particles++, "B0",     3, 5.2794, +0., 1.532e-12);     // G4BMesonZero
  fPDGCode[fNPDGCodes++]=511;          //58 = B0

  Gspart(fNG3Particles++, "B0 bar", 3, 5.2794, -0., 1.532e-12);     // G4AntiBMesonZero
  fPDGCode[fNPDGCodes++]=-511;         //58 = B0bar

  Gspart(fNG3Particles++, "B+",     4, 5.2790, +1., 1.638e-12);     // G4BMesonPlus *
  fPDGCode[fNPDGCodes++]=521;          //60 = B+

  Gspart(fNG3Particles++, "B-",     4, 5.2790, -1., 1.638e-12);     // G4BMesonMinus *
  fPDGCode[fNPDGCodes++]=-521;         //61 = B-

  Gspart(fNG3Particles++, "Bs",     3, 5.3675, +0., 1.466e-12);     // G4BsMesonZero
  fPDGCode[fNPDGCodes++]=531;          //62 = B_s

  Gspart(fNG3Particles++, "Bs bar", 3, 5.3675, -0., 1.466e-12);     // G4AntiBsMesonZero
  fPDGCode[fNPDGCodes++]=-531;         //63 = B_s bar

  Gspart(fNG3Particles++, "Lambda_b",     3, 5.624, +0., 1.24e-12);
  fPDGCode[fNPDGCodes++]=5122;         //64 = Lambda_b

  Gspart(fNG3Particles++, "Lambda_b bar", 3, 5.624, -0., 1.24e-12);
  fPDGCode[fNPDGCodes++]=-5122;        //65 = Lambda_b bar

  Gspart(fNG3Particles++, "J/Psi",       3, 3.096916, 0., 7.6e-21);   // G4JPsi
  fPDGCode[fNPDGCodes++]=443;          // 66 = J/Psi

  Gspart(fNG3Particles++, "Psi Prime",   3, 3.686,   0., 0.);
  fPDGCode[fNPDGCodes++]=20443;        // 67 = Psi prime

  Gspart(fNG3Particles++, "Upsilon(1S)", 3, 9.46037, 0., 0.);
  fPDGCode[fNPDGCodes++]=553;          // 68 = Upsilon(1S)

  Gspart(fNG3Particles++, "Upsilon(2S)", 3, 10.0233, 0., 0.);
  fPDGCode[fNPDGCodes++]=20553;        // 69 = Upsilon(2S)

  Gspart(fNG3Particles++, "Upsilon(3S)", 3, 10.3553, 0., 0.);
  fPDGCode[fNPDGCodes++]=30553;        // 70 = Upsilon(3S)

  Gspart(fNG3Particles++, "Anti Neutrino (e)",       3, 0., 0., 1.e20);
  fPDGCode[fNPDGCodes++]=-12;          // 71 = anti electron neutrino

  Gspart(fNG3Particles++, "Neutrino (mu)",           3, 0., 0., 1.e20);
  fPDGCode[fNPDGCodes++]=14;           // 72 = muon neutrino

  Gspart(fNG3Particles++, "Anti Neutrino (mu)", 3, 0., 0., 1.e20);
  fPDGCode[fNPDGCodes++]=-14;          // 73 = anti muon neutrino

  Gspart(fNG3Particles++, "Neutrino (tau)",     3, 0., 0., 1.e20);
  fPDGCode[fNPDGCodes++]=16;           // 74 = tau neutrino

  Gspart(fNG3Particles++, "Anti Neutrino (tau)",3, 0., 0., 1.e20);
  fPDGCode[fNPDGCodes++]=-16;          // 75 = anti tau neutrino

/* --- Define additional decay modes --- */
/* --- omega(783) --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 33;
    bratio[0] = 89.;
    bratio[1] = 8.5;
    bratio[2] = 2.5;
    mode[0] = 70809;
    mode[1] = 107;
    mode[2] = 908;
    Gsdk(ipa, bratio, mode);
/* --- phi(1020) --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 34;
    bratio[0] = 49.;
    bratio[1] = 34.4;
    bratio[2] = 12.9;
    bratio[3] = 2.4;
    bratio[4] = 1.3;
    mode[0] = 1112;
    mode[1] = 1610;
    mode[2] = 4407;
    mode[3] = 90807;
    mode[4] = 1701;
    Gsdk(ipa, bratio, mode);
/* --- D+ --- */
    /*
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 35;
    bratio[0] = 25.;
    bratio[1] = 25.;
    bratio[2] = 25.;
    bratio[3] = 25.;
    mode[0] = 80809;
    mode[1] = 120808;
    mode[2] = 111208;
    mode[3] = 110809;
    Gsdk(ipa, bratio, mode);
    */
/* --- D- --- */
    /*
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 36;
    bratio[0] = 25.;
    bratio[1] = 25.;
    bratio[2] = 25.;
    bratio[3] = 25.;
    mode[0] = 90908;
    mode[1] = 110909;
    mode[2] = 121109;
    mode[3] = 120908;
    Gsdk(ipa, bratio, mode);
    */
/* --- D0 --- */
    /*
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 37;
    bratio[0] = 33.;
    bratio[1] = 33.;
    bratio[2] = 33.;
    mode[0] = 809;
    mode[1] = 1208;
    mode[2] = 1112;
    Gsdk(ipa, bratio, mode);
    */
/* --- Anti D0 --- */
    /*
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 38;
    bratio[0] = 33.;
    bratio[1] = 33.;
    bratio[2] = 33.;
    mode[0] = 809;
    mode[1] = 1109;
    mode[2] = 1112;
    Gsdk(ipa, bratio, mode);
    */
/* --- rho+ --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 42;
    bratio[0] = 100.;
    mode[0] = 807;
    Gsdk(ipa, bratio, mode);
/* --- rho- --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 43;
    bratio[0] = 100.;
    mode[0] = 907;
    Gsdk(ipa, bratio, mode);
/* --- rho0 --- */
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 44;
    bratio[0] = 100.;
    mode[0] = 809;
    Gsdk(ipa, bratio, mode);
    /*
// --- jpsi ---
    for (kz = 0; kz < 6; ++kz) {
	bratio[kz] = 0.;
	mode[kz] = 0;
    }
    ipa = 113;
    bratio[0] = 50.;
    bratio[1] = 50.;
    mode[0] = 506;
    mode[1] = 605;
    Gsdk(ipa, bratio, mode);
// --- upsilon ---
    ipa = 114;
    Gsdk(ipa, bratio, mode);
// --- phi ---
    ipa = 115;
    Gsdk(ipa, bratio, mode);
    */
//
    AddParticlesToPdgDataBase();
}

//______________________________________________________________________
Int_t TGeant3::VolId(const Text_t *name) const
{
  //
  // Return the unique numeric identifier for volume name
  //
  Int_t gname,i;
  strncpy((char *) &gname, name, 4);
  for(i=1; i<=fGcnum->nvolum; i++)
    if(gname == fZiq[fGclink->jvolum+i]) return i;
  printf("VolId: Volume %s not found\n",name);
  return 0;
}

//______________________________________________________________________
Int_t TGeant3::MediumId(const Text_t *medName) const
{
    // Return the unique numeric identifier for medium name                  

  Int_t nmed = fMedNames.GetEntriesFast();
  for ( Int_t imed = 1; imed < nmed; imed++ ) {
  
    TString name = ((TObjString*)fMedNames.At(imed))->GetString();
    if ( name == TString(medName) )  return imed;
  }
  printf("MediumId: Medium %s not found\n", medName);
  return 0;
}      
        
//______________________________________________________________________
Int_t TGeant3::NofVolumes() const
{
  //
  // Return total number of volumes in the geometry
  //
  return fGcnum->nvolum;
}

//______________________________________________________________________
Int_t TGeant3::NofVolDaughters(const char* volName) const
{
// Return number of daughters of the volume specified by volName
// According to A. Morsch' G3toRoot class
// ---

  Int_t idvol = VolId(volName);

  Int_t jvo = fZlq[fGclink->jvolum-idvol];
  Int_t nin = Int_t(fZq[jvo+3]);
  return nin;
}

//______________________________________________________________________
const char*  TGeant3::VolDaughterName(const char* volName, Int_t i) const
{
// Return the name of i-th daughters of the volume specified by volName
// According to A. Morsch' G3toRoot class
// ---

  Int_t idvol = VolId(volName);

  Int_t jvo = fZlq[fGclink->jvolum-idvol];
  Int_t nin=i+1;
  Int_t jin = fZlq[jvo-nin];
  Int_t idvold = Int_t(fZq[jin+2]);;

  return VolName(idvold);
}


//______________________________________________________________________
Int_t TGeant3::VolDaughterCopyNo(const char* volName, Int_t i) const
{
// Return the copyNo of i-th daughters of the volume specified by volName
// According to A. Morsch' G3toRoot class
// ---

  Int_t idvol = VolId(volName);

  Int_t jvo = fZlq[fGclink->jvolum-idvol];
  Int_t nin=i+1;
  Int_t jin = fZlq[jvo-nin];

  return  Int_t(fZq[jin +3]);
}

//______________________________________________________________________
Int_t TGeant3::VolId2Mate(Int_t id) const
{
  //
  // Return material number for a given volume id
  //
  if(id<1 || id > fGcnum->nvolum || fGclink->jvolum<=0)
    return 0;
  else {
    Int_t jvo = fZlq[fGclink->jvolum-id];
    return Int_t(fZq[jvo+4]);
  }
}

//______________________________________________________________________
const char* TGeant3::VolName(Int_t id) const
{
  //
  // Return the volume name given the volume identifier
  //
  if(id<1 || id > fGcnum->nvolum || fGclink->jvolum<=0)
    return fVolNames[fGcnum->nvolum];
  else
    return fVolNames[id-1];
}

//______________________________________________________________________
Bool_t  TGeant3::SetCut(const char* cutName, Double_t cutValue)
{
  //
  // Set transport cuts for particles
  //
  Bool_t success = kTRUE;

  if(!strcmp(cutName,"CUTGAM"))
    fGccuts->cutgam=cutValue;
  else if(!strcmp(cutName,"CUTELE"))
    fGccuts->cutele=cutValue;
  else if(!strcmp(cutName,"CUTNEU"))
    fGccuts->cutneu=cutValue;
  else if(!strcmp(cutName,"CUTHAD"))
    fGccuts->cuthad=cutValue;
  else if(!strcmp(cutName,"CUTMUO"))
    fGccuts->cutmuo=cutValue;
  else if(!strcmp(cutName,"BCUTE"))
    fGccuts->bcute=cutValue;
  else if(!strcmp(cutName,"BCUTM"))
    fGccuts->bcutm=cutValue;
  else if(!strcmp(cutName,"DCUTE"))
    fGccuts->dcute=cutValue;
  else if(!strcmp(cutName,"DCUTM"))
    fGccuts->dcutm=cutValue;
  else if(!strcmp(cutName,"PPCUTM"))
    fGccuts->ppcutm=cutValue;
  else if(!strcmp(cutName,"TOFMAX"))
    fGccuts->tofmax=cutValue;
  else {
    Warning("SetCut","Cut %s not implemented\n",cutName);
    success = kFALSE;
  }

  return success;
}

//______________________________________________________________________
Bool_t  TGeant3::SetProcess(const char* flagName, Int_t flagValue)
{
  //
  // Set thresholds for different processes
  //
  Bool_t success = kTRUE;

  if(!strcmp(flagName,"PAIR"))
    fGcphys->ipair=flagValue;
  else if(!strcmp(flagName,"COMP"))
    fGcphys->icomp=flagValue;
  else if(!strcmp(flagName,"PHOT"))
    fGcphys->iphot=flagValue;
  else if(!strcmp(flagName,"PFIS"))
    fGcphys->ipfis=flagValue;
  else if(!strcmp(flagName,"DRAY"))
    fGcphys->idray=flagValue;
  else if(!strcmp(flagName,"ANNI"))
    fGcphys->ianni=flagValue;
  else if(!strcmp(flagName,"BREM"))
    fGcphys->ibrem=flagValue;
  else if(!strcmp(flagName,"HADR"))
    fGcphys->ihadr=flagValue;
  else if(!strcmp(flagName,"MUNU"))
    fGcphys->imunu=flagValue;
  else if(!strcmp(flagName,"DCAY"))
    fGcphys->idcay=flagValue;
  else if(!strcmp(flagName,"LOSS"))
    fGcphys->iloss=flagValue;
  else if(!strcmp(flagName,"MULS"))
    fGcphys->imuls=flagValue;
  else if(!strcmp(flagName,"RAYL"))
    fGcphys->irayl=flagValue;
  else if(!strcmp(flagName,"STRA"))
    fGcphlt->istra=flagValue;
  else if(!strcmp(flagName,"SYNC"))
    fGcphlt->isync=flagValue;
  else if(!strcmp(flagName,"CKOV"))
    fGctlit->itckov = flagValue;
  else  {
    Warning("SetFlag","Flag %s not implemented\n",flagName);
    success = kFALSE;
  }

  return  success;
}

 //______________________________________________________________________
Bool_t TGeant3::DefineParticle(Int_t pdg,const char* name,TMCParticleType type,
                      Double_t mass, Double_t charge, Double_t lifetime)
{
// Old function definition, now replaced with more arguments

  TVirtualMC::DefineParticle(pdg, name, type, mass, charge, lifetime);
  
  return false;
}                        
                      

//______________________________________________________________________
Bool_t TGeant3::DefineParticle(Int_t pdg,const char* name, TMCParticleType mcType,
                      Double_t mass, Double_t charge, Double_t lifetime,
                      const TString& /*pType*/, Double_t /*width*/, 
                      Int_t /*iSpin*/, Int_t /*iParity*/, Int_t /*iConjugation*/, 
                      Int_t /*iIsospin*/, Int_t /*iIsospinZ*/, Int_t /*gParity*/,
                      Int_t /*lepton*/, Int_t /*baryon*/,
                      Bool_t /*stable*/, Bool_t /*shortlived*/,
                      const TString& /*subType*/,
                      Int_t /*antiEncoding*/, Double_t /*magMoment*/,
                      Double_t /*excitation*/)
{
//
// Set a user defined particle
// Function is ignored if particle with specified pdg
// already exists and error report is printed.
// ---

  // Check if particle with specified pdg already exists
  // in TGeant3
  if (IdFromPDG(pdg) > 0) {
    Error("SetParticle", "Particle already exists.");
    return kFALSE;
  }

  // Check if particle type is known to Geant3
  Int_t itrtyp = TransportMethod(mcType);
  if (itrtyp < 0) {
    Error("SetParticle", "Unknown particle transport.");
    return kFALSE;
  }

  // Add particle to Geant3
  Gspart(fNG3Particles++, name, itrtyp, mass, charge, lifetime);

  // Add particle to TDatabasePDG
  // (if it does not yet exist here)
  if (!TDatabasePDG::Instance()->GetParticle(pdg))
    TDatabasePDG::Instance()
      ->AddParticle(name, name, mass, kTRUE, 0, charge*3,
                    ParticleClass(mcType).Data(), pdg);
                    
  // Resize fPDGCode table if needed
  if ( fNPDGCodes >= fPDGCode.GetSize() ) 
    fPDGCode.Set( fPDGCode.GetSize() + 100);                 

  fPDGCode[fNPDGCodes++] = pdg;

  return kTRUE;
}

//______________________________________________________________________
Bool_t  TGeant3::DefineIon(const char* name, Int_t Z, Int_t A, Int_t Q,
                           Double_t /* excEnergy */, Double_t mass)
{
//
// Set a user defined ion.
// ---

  // Define pdgEncoding
  //
  Int_t pdg = GetIonPdg(Z, A);
  Int_t pdgMax = pdg + 9;

  // Find isomer number which is not yet used
  while (TDatabasePDG::Instance()->GetParticle(pdg) &&
         pdg < pdgMax)
      pdg++;
  if (TDatabasePDG::Instance()->GetParticle(pdg)) {
      Fatal("SetIon", "All isomer numbers are already used");
      return kFALSE;
  }

  // Particle properties
        // excitation energy not used by G3
  if (mass < 1e-09) mass = 0.9382723 * A;
     // approximative mass if not specified by user
  Double_t charge = Q;
  TMCParticleType partType = kPTIon;
  Double_t lifetime = 1.e20;

  // Call DefineParticle now
  return DefineParticle(
           pdg, name, partType, mass, charge, lifetime,
           "nucleus", 0.0, 1, 1, 0, 1, 1, 0, 0, 1, kTRUE);
}

//______________________________________________________________________
TString  TGeant3::ParticleName(Int_t pdg) const
{
//  Return G3 particle name
// ---

  char name[21];
  Int_t itrtyp;
  Float_t amass, charge, tlife;
  Gfpart(pdg, name, itrtyp,amass, charge, tlife);
  name[20] = '\0';

  return TString(name);
}

//______________________________________________________________________
Double_t  TGeant3::ParticleMass(Int_t pdg) const
{
//  Return G3 particle mass
// ---

  char name[20];
  Int_t itrtyp;
  Float_t mass, charge, tlife;
  Gfpart(pdg,name, itrtyp, mass, charge, tlife);

  return mass;
}

//______________________________________________________________________
Double_t  TGeant3::ParticleCharge(Int_t pdg) const
{
// Return G3 particle charge (in e)
// ---

  char name[20];
  Int_t itrtyp;
  Float_t mass, charge, tlife;
  Gfpart(pdg,name, itrtyp, mass, charge, tlife);

  return charge;
}

//______________________________________________________________________
Double_t  TGeant3::ParticleLifeTime(Int_t pdg) const
{
// Return G3 particle life time
// ---

  char name[20];
  Int_t itrtyp;
  Float_t mass, charge, tlife;
  Gfpart(pdg, name, itrtyp, mass, charge, tlife);

  return tlife;
}

//______________________________________________________________________
TMCParticleType TGeant3::ParticleMCType(Int_t pdg) const
{
// Return MC particle type
// ---

  char name[20];
  Int_t itrtyp;
  Float_t mass, charge, tlife;
  Gfpart(pdg,name, itrtyp, mass, charge, tlife);

  return ParticleType(itrtyp);
}


//______________________________________________________________________
Double_t TGeant3::Xsec(char* reac, Double_t /* energy */,
		      Int_t part, Int_t /* mate */)
{
  //
  // Calculate X-sections -- dummy for the moment
  //
  if(!strcmp(reac,"PHOT"))
  {
    if(part!=22) {
      Error("Xsec","Can calculate photoelectric only for photons\n");
    }
  }
  return 0;
}

//______________________________________________________________________
void TGeant3::TrackPosition(TLorentzVector &xyz) const
{
  //
  // Return the current position in the master reference frame of the
  // track being transported
  //
  xyz[0]=fGctrak->vect[0];
  xyz[1]=fGctrak->vect[1];
  xyz[2]=fGctrak->vect[2];
  xyz[3]=fGctrak->tofg;
}

//______________________________________________________________________
void TGeant3::TrackPosition(Double_t &x, Double_t &y, Double_t &z) const
{
  //
  // Return the current position in the master reference frame of the
  // track being transported
  //
  x=fGctrak->vect[0];
  y=fGctrak->vect[1];
  z=fGctrak->vect[2];
}

//______________________________________________________________________
Double_t TGeant3::TrackTime() const
{
  //
  // Return the current time of flight of the track being transported
  //
  return fGctrak->tofg;
}

//______________________________________________________________________
void TGeant3::TrackMomentum(TLorentzVector &xyz) const
{
  //
  // Return the direction and the momentum (GeV/c) of the track
  // currently being transported
  //
  Double_t ptot=fGctrak->vect[6];
  xyz[0]=fGctrak->vect[3]*ptot;
  xyz[1]=fGctrak->vect[4]*ptot;
  xyz[2]=fGctrak->vect[5]*ptot;
  xyz[3]=fGctrak->getot;
}

//______________________________________________________________________
void TGeant3::TrackMomentum(Double_t &px, Double_t &py, Double_t &pz, 
                            Double_t &etot) const
{
  //
  // Return the direction and the momentum (GeV/c) of the track
  // currently being transported
  //
  Double_t ptot=fGctrak->vect[6];
  px  =fGctrak->vect[3]*ptot;
  py  =fGctrak->vect[4]*ptot;
  pz  =fGctrak->vect[5]*ptot;
  etot=fGctrak->getot;
}

//______________________________________________________________________
Double_t TGeant3::TrackCharge() const
{
  //
  // Return charge of the track currently transported
  //
  return fGckine->charge;
}

//______________________________________________________________________
Double_t TGeant3::TrackMass() const
{
  //
  // Return the mass of the track currently transported
  //
  return fGckine->amass;
}

//______________________________________________________________________
Int_t TGeant3::TrackPid() const
{
  //
  // Return the id of the particle transported
  //
  return PDGFromId(fGckine->ipart);
}

//______________________________________________________________________
Double_t TGeant3::TrackStep() const
{
  //
  // Return the length in centimeters of the current step
  //
  return fGctrak->step;
}

//______________________________________________________________________
Double_t TGeant3::TrackLength() const
{
  //
  // Return the length of the current track from its origin
  //
  return fGctrak->sleng;
}

//______________________________________________________________________
Bool_t TGeant3::IsNewTrack() const
{
  //
  // True if the track is not at the boundary of the current volume
  //
  return (fGctrak->sleng==0);
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackInside() const
{
  //
  // True if the track is not at the boundary of the current volume
  //
  return (fGctrak->inwvol==0);
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackEntering() const
{
  //
  // True if this is the first step of the track in the current volume
  //
  return (fGctrak->inwvol==1);
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackExiting() const
{
  //
  // True if this is the last step of the track in the current volume
  //
  return (fGctrak->inwvol==2);
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackOut() const
{
  //
  // True if the track is out of the setup
  //
  return (fGctrak->inwvol==3);
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackStop() const
{
  //
  // True if the track energy has fallen below the threshold
  //
  return (fGctrak->istop==2);
}

//______________________________________________________________________
Int_t   TGeant3::NSecondaries() const
{
  //
  // Number of secondary particles generated in the current step
  //
  return fGcking->ngkine;
}

//______________________________________________________________________
Int_t   TGeant3::CurrentEvent() const
{
  //
  // Number of the current event
  //
  return fGcflag->idevt;
}

//______________________________________________________________________
TMCProcess TGeant3::ProdProcess(Int_t ) const
{
  //
  // Name of the process that has produced the secondary particles
  // in the current step

  //  Modified: to make use of GCKING/KCASE variable for determining the production
  //  mechanism of the secondaries.  The old method was to pick the first 
  //  active process from the current step's list of active processes 
  //  that had the capability of generating secondaries.  This occasionally 
  //  picked the wrong secondary production mechanism.

  Int_t imech=0;

  if ( fGcking->ngkine <= 0 ) return kPNoProcess;

  // Secondaries generated, determine production mechanism hollerith 
  for (Int_t km = 0; km < MAXMEC; ++km) {
    if ( fGcking->kcase == fGctrak->namec[km] ) {
      imech = km;
      break;
    }
  }

  TMCProcess vmcmech = G3toVMC(imech+1);
  if ( vmcmech == kPNoProcess ) {
    // failure to find matching process
    printf(
    "* TGeant3::ProdProcess secondaries present,but no matching process!* \n");
  }

  return vmcmech;
}

//______________________________________________________________________
Int_t TGeant3::StepProcesses(TArrayI &proc) const
{
  //
  // Return processes active in the current step
  //
  Int_t i;
  Int_t nproc=Gctrak()->nmec;
  
  // Set no active process if there are no processes
  if (nproc==0) {
    proc.Set(1);
    proc[0] = kPNull;
    return 1;
  }  
  
  //
  proc.Set(nproc);
  Int_t nvproc=0;
  //
  for (i=0; i<nproc; ++i)
    if((proc[nvproc]=G3toVMC(Gctrak()->lmec[i]))!=kPNoProcess) nvproc++;
  //
  proc.Set(nvproc);
  //
  return nvproc;
}

//______________________________________________________________________
TMCProcess TGeant3::G3toVMC(Int_t iproc) const
{
  //
  // Conversion between GEANT and TMC processes
  //

  const TMCProcess kPG2MC1[30] = {
    kPTransportation, kPMultipleScattering, kPEnergyLoss, kPMagneticFieldL, kPDecay,
    kPPair, kPCompton, kPPhotoelectric, kPBrem, kPDeltaRay,
    kPAnnihilation, kPHadronic, kPHCElastic, kPEvaporation, kPNuclearFission,
    kPNuclearAbsorption, kPPbarAnnihilation, kPNCapture, kPHIElastic, 
    kPHInhelastic, kPMuonNuclear, kPTOFlimit, kPPhotoFission, kPNoProcess, 
    kPRayleigh, kPNoProcess, kPNoProcess, kPNoProcess, kPNull, kPStop};

  const TMCProcess kPG2MC2[9] = {
      kPLightAbsorption, kPLightScattering, kStepMax, kPNoProcess, kPCerenkov,
      kPLightReflection, kPLightRefraction, kPSynchrotron, kPNoProcess};

  TMCProcess proc=kPNoProcess;
  if(0<iproc && iproc<=30) proc= kPG2MC1[iproc-1];
  else if(101<=iproc && iproc<=109) proc= kPG2MC2[iproc-100-1];
  return proc;
}


//______________________________________________________________________
void    TGeant3::GetSecondary(Int_t isec, Int_t& ipart,
			      TLorentzVector &x, TLorentzVector &p)
{
  //
  // Get the parameters of the secondary track number isec produced
  // in the current step
  //
  Int_t i;
  if(-1<isec && isec<fGcking->ngkine) {
    ipart=Int_t (fGcking->gkin[isec][4] +0.5);
    for(i=0;i<3;i++) {
      x[i]=fGckin3->gpos[isec][i];
      p[i]=fGcking->gkin[isec][i];
    }
    x[3]=fGcking->tofd[isec];
    p[3]=fGcking->gkin[isec][3];
  } else {
    printf(" * TGeant3::GetSecondary * Secondary %d does not exist\n",isec);
    x[0]=x[1]=x[2]=x[3]=p[0]=p[1]=p[2]=p[3]=0;
    ipart=0;
  }
}

//______________________________________________________________________
void TGeant3::InitLego()
{
  //
  // Set switches for lego transport
  //
  SetSWIT(4,0);
  SetDEBU(0,0,0);  //do not print a message
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackDisappeared() const
{
  //
  // True if the current particle has disappeared
  // either because it decayed or because it underwent
  // an inelastic collision
  //
  return (fGctrak->istop==1);
}

//______________________________________________________________________
Bool_t TGeant3::IsTrackAlive() const
{
  //
  // True if the current particle is alive and will continue to be
  // transported
  //
  return (fGctrak->istop==0);
}

//______________________________________________________________________
void TGeant3::StopTrack()
{
  //
  // Stop the transport of the current particle and skip to the next
  //
  fGctrak->istop=1;
}

//______________________________________________________________________
void TGeant3::StopEvent()
{
  //
  // Stop simulation of the current event and skip to the next
  //
  fGcflag->ieotri=1;
}

//______________________________________________________________________
void TGeant3::StopRun()
{
  //
  // Stop simulation of the current event and set the abort run flag to true
  //

  StopTrack();
  StopEvent();
  fStopRun = kTRUE;
}

//______________________________________________________________________
Double_t TGeant3::MaxStep() const
{
  //
  // Return the maximum step length in the current medium
  //
  return fGctmed->stemax;
}

//______________________________________________________________________
void TGeant3::SetMaxStep(Double_t maxstep)
{
  //
  // Set the maximum step allowed till the particle is in the current medium
  //
  fGctmed->stemax=maxstep;
}

//______________________________________________________________________
void TGeant3::SetMaxNStep(Int_t maxnstp)
{
  //
  // Set the maximum number of steps till the particle is in the current medium
  //
  fGctrak->maxnst=maxnstp;
}

void TGeant3::ForceDecayTime(Float_t time)
{
    //
    // Force the decay time of the current particle
    //
    TLorentzVector p;
    TrackMomentum(p);
    Gcphys()->sumlif = time / p.Beta() / p.Gamma()  * 2.99792458e10;
}

//______________________________________________________________________
Int_t TGeant3::GetMaxNStep() const
{
  //
  // Maximum number of steps allowed in current medium
  //
  return fGctrak->maxnst;
}

//_______________________________________________________________________
void TGeant3::G3Material(Int_t& kmat, const char* name, Double_t a, 
                         Double_t z, Double_t dens, Double_t radl, 
                         Double_t absl, Float_t* buf, Int_t nwbuf)
{
  //
  // Defines a Material
  //
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorption length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //
  Int_t jmate=fGclink->jmate;
  kmat=1;
  Int_t ns, i;
  if(jmate>0) {
    ns=fZiq[jmate-2];
    kmat=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jmate-i]==0) {
	kmat=i;
	break;
      }
    }
  }
  Float_t fa = a;
  Float_t fz = z;
  Float_t fdens = dens;
  Float_t fradl = radl;
  Float_t fabsl = absl;

  g3smate(kmat,PASSCHARD(name), fa,  fz, fdens, fradl, fabsl, buf,
	 nwbuf PASSCHARL(name));
}

//______________________________________________________________________
void TGeant3::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
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
  //  absl               absorption length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //

  G3Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
}

//______________________________________________________________________
void TGeant3::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
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
  //  absl               absorption length in cm
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
}

//______________________________________________________________________
void TGeant3::G3Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z,
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
  // weights.
  //

  Int_t jmate=fGclink->jmate;
  kmat=1;
  Int_t ns, i;
  if(jmate>0) {
    ns=fZiq[jmate-2];
    kmat=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jmate-i]==0) {
	kmat=i;
	break;
      }
    }
  }
  g3smixt(kmat,PASSCHARD(name),a,z,Float_t(dens),nlmat,wmat PASSCHARL(name));
}

//______________________________________________________________________
void TGeant3::Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z,
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
  // weights.
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
}

//______________________________________________________________________
void TGeant3::Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z,
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
  // weights.
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
}

//______________________________________________________________________
void TGeant3::G3Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
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
  Int_t jtmed=fGclink->jtmed;
  kmed=1;
  Int_t ns, i;
  if(jtmed>0) {
    ns=fZiq[jtmed-2];
    kmed=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jtmed-i]==0) {
	kmed=i;
	break;
      }
    }
  }
  Float_t ffieldm = fieldm;
  Float_t ftmaxfd = tmaxfd;
  Float_t fstemax = stemax;
  Float_t fdeemax = deemax;
  Float_t fepsil  = epsil;
  Float_t fstmin =  stmin;
  g3stmed(kmed, PASSCHARD(name),nmat,isvol,ifield,ffieldm,ftmaxfd,fstemax,
          fdeemax, fepsil, fstmin, ubuf, nbuf PASSCHARL(name));

  fMedNames.AddAtAndExpand(new TObjString(name), kmed);           
}

//______________________________________________________________________
void TGeant3::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
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

  G3Medium(kmed,name,nmat,isvol,ifield,fieldm,tmaxfd,stemax,deemax,epsil,
           stmin, ubuf, nbuf);

}

//______________________________________________________________________
void TGeant3::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
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
  //  stmin     min. step due to continuous processes (cm)
  //
  //  ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
  //  ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
  //  performed with g3helix; ifield = 3 if tracking performed with g3helx3.
  //

  Float_t* fubuf = CreateFloatArray(ubuf, nbuf);
  G3Medium(kmed,name,nmat,isvol,ifield,fieldm,tmaxfd,stemax,deemax,epsil,
           stmin, fubuf, nbuf);
  delete [] fubuf;

}

//______________________________________________________________________
void TGeant3::Matrix(Int_t& krot, Double_t thex, Double_t phix, Double_t they,
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
  Int_t jrotm=fGclink->jrotm;
  krot=1;
  Int_t ns, i;
  if(jrotm>0) {
    ns=fZiq[jrotm-2];
    krot=ns+1;
    for(i=1; i<=ns; i++) {
      if(fZlq[jrotm-i]==0) {
	krot=i;
	break;
      }
    }
  }
  g3srotm(krot, thex, phix, they, phiy, thez, phiz);
}

//______________________________________________________________________
Int_t TGeant3::CurrentMedium() const
{
  //
  // Return the number of the current medium
  //
//#ifdef WITHROOT
//  Int_t imed = 0;
//  TGeoNode *node = gGeoManager->GetCurrentNode();
//  if (!node) imed = gGeoManager->GetTopNode()->GetVolume()->
//                                                 GetMedium()->GetId();
//  else       imed = node->GetVolume()->GetMedium()->GetId();
  //printf("==GetMedium: ROOT id=%i  numed=%i\n", imed,fGctmed->numed);
//#endif
  return fGctmed->numed;
}

//______________________________________________________________________
Int_t TGeant3::GetMedium() const
{
  //
  // Return the number of the current medium
  // Deprecated function - replaced with CurrentMedium()
  //

  Warning("GetMedium", 
          "Deprecated function - use CurrentMedium() instead");

  return CurrentMedium();
}
 
//______________________________________________________________________
void  TGeant3::SetRootGeometry()
{
// Notify Geant3 about use of TGeo geometry.
// The materials and tracking medias will be imported from
// TGeo at FinishGeometry().

  Fatal("SetRootGeometry",
        "TGeant3 does not support Root geometry");

  fImportRootGeometry = kTRUE;
}  

//______________________________________________________________________
void TGeant3::SetUserParameters(Bool_t isUserParameters)
{
// Activate the parameters defined in tracking media
// (DEEMAX, STMIN, STEMAX), which are, be default, ignored.

  SetAUTO(!isUserParameters);
}  

//______________________________________________________________________
const char *TGeant3::GetPath()
{
// Get current path inside G3 geometry
   Int_t i,j;
   if ((i=fGcvolu->nlevel-1)<0) {
      Warning("GetPath", "level null");
      return fPath;
   }
   fPath[0] = '/';
   char name[10];
   char *namcur = fPath+1;
   Int_t gname, copy;
   Int_t nch=0;
   for (j=0; j<i+1; j++) {
      gname = fGcvolu->names[j];
      copy = fGcvolu->number[j];
      memcpy(name, &gname, 4);
      name[4]=0;
      sprintf(namcur, "%s_%d/", name, copy);
      nch = strlen(fPath);
      namcur = fPath+nch;
   }
   fPath[nch-1]=0;
   return fPath;
}

//______________________________________________________________________
const char *TGeant3::GetNodeName()
{
// Get name of current G3 node
   Int_t i=fGcvolu->nlevel-1;
   if (i<0) return "";
   Int_t gname = fGcvolu->names[i];
   Int_t copy = fGcvolu->number[i];
   char name[10];
   memcpy(name, &gname, 4);
   name[4] = 0;
   sprintf(fPath, "%s_%d", name, copy);
   return fPath;
}

//______________________________________________________________________
Double_t TGeant3::Edep() const
{
  //
  // Return the energy lost in the current step
  //
  return fGctrak->destep;
}

//______________________________________________________________________
Double_t TGeant3::Etot() const
{
  //
  // Return the total energy of the current track
  //
  return fGctrak->getot;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GBASE
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gfile(const char * /*filename*/, const char * /*option*/)
{
  //
  //    Routine to open a GEANT/RZ data base.
  //
  //    LUN logical unit number associated to the file
  //
  //    CHFILE RZ file name
  //
  //    CHOPT is a character string which may be
  //        N  To create a new file
  //        U  to open an existing file for update
  //       " " to open an existing file for read only
  //        Q  The initial allocation (default 1000 records)
  //           is given in IQUEST(10)
  //        X  Open the file in exchange format
  //        I  Read all data structures from file to memory
  //        O  Write all data structures from memory to file
  //
  // Note:
  //      If options "I"  or "O" all data structures are read or
  //         written from/to file and the file is closed.
  //      See routine GRMDIR to create subdirectories
  //      See routines GROUT,GRIN to write,read objects
  //
  //g3rfile(21, PASSCHARD(filename), PASSCHARD(option) PASSCHARL(filename)
//	 PASSCHARL(option));
}

//______________________________________________________________________
void  TGeant3::Gpcxyz()
{
  //
  //    Print track and volume parameters at current point
  //

    g3pcxyz();
}
//______________________________________________________________________
void  TGeant3::Ggclos()
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
  g3gclos();
  // Create internal list of volumes
  fVolNames = new char[fGcnum->nvolum+1][5];
  Int_t i;
  for(i=0; i<fGcnum->nvolum; ++i) {
    strncpy(fVolNames[i], (char *) &fZiq[fGclink->jvolum+i+1], 4);
    fVolNames[i][4]='\0';
  }
  strcpy(fVolNames[fGcnum->nvolum],"NULL");
}

//______________________________________________________________________
void  TGeant3::Glast()
{
  //
  // Finish a Geant run
  //
  g3last();
}

//______________________________________________________________________
void  TGeant3::Gprint(const char *name)
{
  //
  // Routine to print data structures
  // CHNAME   name of a data structure
  //
  char vname[5];
  Vname(name,vname);
  g3print(PASSCHARD(vname),0 PASSCHARL(vname));
}

//______________________________________________________________________
void  TGeant3::Grun()
{
  //
  // Steering function to process one run
  //
  g3run();
}

//______________________________________________________________________
void  TGeant3::Gtrig()
{
  //
  // Steering function to process one event
  //
  g3trig();

  //printf("count_gmedia= %8d\n",count_gmedia);
  //printf("count_gtmedi= %8d\n",count_gtmedi);
  //printf("count_ginvol= %8d\n",count_ginvol);
  //printf("count_gtnext= %8d\n",count_gtnext);
}

//______________________________________________________________________
void  TGeant3::Gtrigc()
{
  //
  // Clear event partition
  //
  g3trigc();
}

//______________________________________________________________________
void  TGeant3::Gtrigi()
{
  //
  // Initializes event partition
  //
  g3trigi();
}

//______________________________________________________________________
void  TGeant3::Gwork(Int_t nwork)
{
  //
  // Allocates workspace in ZEBRA memory
  //
  g3work(nwork);
}

//______________________________________________________________________
void  TGeant3::Gzinit()
{
  //
  // To initialize GEANT/ZEBRA data structures
  //
  g3zinit();
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,
		      Float_t &dens, Float_t &radl, Float_t &absl,
		      Float_t* ubuf, Int_t& nbuf)
{
  //
  // Return parameters for material IMAT
  //
  g3fmate(imat, PASSCHARD(name), a, z, dens, radl, absl, ubuf, nbuf
	 PASSCHARL(name));
}

//______________________________________________________________________
void  TGeant3::Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,
		      Double_t &dens, Double_t &radl, Double_t &absl,
		      Double_t* ubuf, Int_t& nbuf)
{
  //
  // Return parameters for material IMAT
  //
  Float_t fa = a;
  Float_t fz = z;
  Float_t fdens = dens;
  Float_t fradl = radl;
  Float_t fabsl = absl;
  Float_t* fubuf = CreateFloatArray(ubuf, nbuf);

  Gfmate(imat, name, fa, fz, fdens, fradl, fabsl, fubuf, nbuf);

  a = fa;
  z = fz;
  dens = fdens;
  radl = fradl;
  absl = fabsl;
  for (Int_t i=0; i<nbuf; i++) ubuf[i] = fubuf[i];

  delete [] fubuf;
}

//______________________________________________________________________
void  TGeant3::Gfpart(Int_t ipart, char *name, Int_t &itrtyp,
		   Float_t &amass, Float_t &charge, Float_t &tlife) const
{
  //
  // Return parameters for particle of type IPART
  //
  //Float_t *ubuf=0;
  Float_t ubuf[100];
  Int_t   nbuf;
  Int_t igpart = IdFromPDG(ipart);
  g3fpart(igpart, PASSCHARD(name), itrtyp, amass, charge, tlife, ubuf, nbuf
	 PASSCHARL(name));
}

//______________________________________________________________________
void  TGeant3::Gftmed(Int_t numed, char *name, Int_t &nmat, Int_t &isvol,
		   Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd,
		    Float_t &stemax, Float_t &deemax, Float_t &epsil,
		    Float_t &stmin, Float_t *ubuf, Int_t *nbuf)
{
  //
  // Return parameters for tracking medium NUMED
  //
  g3ftmed(numed, PASSCHARD(name), nmat, isvol, ifield, fieldm, tmaxfd, stemax,
         deemax, epsil, stmin, ubuf, nbuf PASSCHARL(name));
}


//______________________________________________________________________
 void  TGeant3::Gftmat(Int_t imate, Int_t ipart, char *chmeca, Int_t kdim,
		      Float_t* tkin, Float_t* value, Float_t* pcut,
		      Int_t &ixst)
{
  //
  // Return parameters for material imate
  //
  g3ftmat(imate, ipart, PASSCHARD(chmeca), kdim,
	 tkin, value, pcut, ixst PASSCHARL(chmeca));

}

//______________________________________________________________________
Float_t TGeant3::Gbrelm(Float_t z, Float_t t, Float_t bcut)
{
  //
  // To calculate energy loss due to soft muon BREMSSTRAHLUNG
  //
  return g3brelm(z,t,bcut);
}

//______________________________________________________________________
Float_t TGeant3::Gprelm(Float_t z, Float_t t, Float_t bcut)
{
  //
  // To calculate DE/DX in GeV*barn/atom for direct pair production by muons
  //
  return g3prelm(z,t,bcut);
}

//______________________________________________________________________
void  TGeant3::Gmate()
{
  //
  // Define standard GEANT materials
  //
  g3mate();
}

//______________________________________________________________________
void  TGeant3::Gpart()
{
  //
  //  Define standard GEANT particles plus selected decay modes
  //  and branching ratios.
  //
  g3part();
}

//______________________________________________________________________
void  TGeant3::Gsdk(Int_t ipart, Float_t *bratio, Int_t *mode)
{
//  Defines branching ratios and decay modes for standard
//  GEANT particles.
   g3sdk(ipart,bratio,mode);
}

//______________________________________________________________________
void  TGeant3::Gsmate(Int_t imat, const char *name, Float_t a, Float_t z,
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
  //  absl               absorption length in cm
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
     Warning("Gsmate","Density was o, set to 0.01 for imat=%d, name=%s",
             imat,name);
     dens = 0.01;
  }
  g3smate(imat,PASSCHARD(name), a, z, dens, radl, absl, ubuf, nbuf
	 PASSCHARL(name));
}

//______________________________________________________________________
void  TGeant3::Gsmixt(Int_t imat, const char *name, Float_t *a, Float_t *z,
		   Float_t dens, Int_t nlmat, Float_t *wmat)
{
  //
  //       Defines mixture OR COMPOUND IMAT as composed by
  //       THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  //
  //       If NLMAT.GT.0 then WMAT contains the PROPORTION BY
  //       WEIGHTS OF EACH BASIC MATERIAL IN THE MIXTURE.
  //
  //       If NLMAT.LT.0 then WMAT contains the number of atoms
  //       of a given kind into the molecule of the COMPOUND
  //       In this case, WMAT in output is changed to relative
  //       weights.
  //
  g3smixt(imat,PASSCHARD(name), a, z,dens, nlmat,wmat PASSCHARL(name));
}

//______________________________________________________________________
void  TGeant3::Gspart(Int_t ipart, const char *name, Int_t itrtyp,
		   Double_t amass, Double_t charge, Double_t tlife)
{
  //
  // Store particle parameters
  //
  // ipart           particle code
  // name            particle name
  // itrtyp          transport method (see GEANT manual)
  // amass           mass in GeV/c2
  // charge          charge in electron units
  // tlife           lifetime in seconds
  //
  Float_t *ubuf=0;
  Int_t   nbuf=0;
  Float_t fmass = amass;
  Float_t fcharge = charge;
  Float_t flife = tlife;

  g3spart(ipart,PASSCHARD(name), itrtyp, fmass, fcharge, flife, ubuf, nbuf
	 PASSCHARL(name));
}

//______________________________________________________________________
void  TGeant3::Gstmed(Int_t numed, const char *name, Int_t nmat, Int_t isvol,
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
}

//______________________________________________________________________
void  TGeant3::Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			   Float_t *absco, Float_t *effic, Float_t *rindex)
{
  //
  //    Stores the tables for UV photon tracking in medium ITMED
  //    Please note that it is the user's responsibility to
  //    provide all the coefficients:
  //
  //
  //       ITMED       Tracking medium number
  //       NPCKOV      Number of bins of each table
  //       PPCKOV      Value of photon momentum (in GeV)
  //       ABSCO       Absorption coefficients
  //                   dielectric: absorption length in cm
  //                   metals    : absorption fraction (0<=x<=1)
  //       EFFIC       Detection efficiency for UV photons
  //       RINDEX      Refraction index (if=0 metal)
  //
  g3sckov(itmed,npckov,ppckov,absco,effic,rindex);
}

//______________________________________________________________________
void  TGeant3::SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			   Float_t *absco, Float_t *effic, Float_t *rindex)
{
  //
  //    Stores the tables for UV photon tracking in medium ITMED
  //    Please note that it is the user's responsibility to
  //    provide all the coefficients:
  //
  //
  //       ITMED       Tracking medium number
  //       NPCKOV      Number of bins of each table
  //       PPCKOV      Value of photon momentum (in GeV)
  //       ABSCO       Absorption coefficients
  //                   dielectric: absorption length in cm
  //                   metals    : absorption fraction (0<=x<=1)
  //       EFFIC       Detection efficiency for UV photons
  //       RINDEX      Refraction index (if=0 metal)
  //
  g3sckov(itmed,npckov,ppckov,absco,effic,rindex);
}

//______________________________________________________________________
void  TGeant3::SetCerenkov(Int_t itmed, Int_t npckov, Double_t *ppckov,
			   Double_t *absco, Double_t *effic, Double_t *rindex)
{
  //
  //    Stores the tables for UV photon tracking in medium ITMED
  //    Please note that it is the user's responsibility to
  //    provide all the coefficients:
  //
  //
  //       ITMED       Tracking medium number
  //       NPCKOV      Number of bins of each table
  //       PPCKOV      Value of photon momentum (in GeV)
  //       ABSCO       Absorption coefficients
  //                   dielectric: absorption length in cm
  //                   metals    : absorption fraction (0<=x<=1)
  //       EFFIC       Detection efficiency for UV photons
  //       RINDEX      Refraction index (if=0 metal)
  //

  Float_t* fppckov = CreateFloatArray(ppckov, npckov);
  Float_t* fabsco  = CreateFloatArray(absco, npckov);
  Float_t* feffic  = CreateFloatArray(effic, npckov);
  Float_t* frindex = CreateFloatArray(rindex, npckov);

  SetCerenkov(itmed, npckov, fppckov, fabsco, feffic, frindex);

  delete [] fppckov;
  delete [] fabsco;
  delete [] feffic;
  delete [] frindex;
}

//______________________________________________________________________
void  TGeant3::DefineOpSurface(const char* name,
                EMCOpSurfaceModel /*model*/, EMCOpSurfaceType /*surfaceType*/,
                EMCOpSurfaceFinish /*surfaceFinish*/, Double_t /*sigmaAlpha*/) 
{
   
   Warning("DefineOpSurface", 
           Form("Called for surface %s. Not applicable in Geant3 - setting is ignored.", name));
}   
                
//______________________________________________________________________
void  TGeant3::SetBorderSurface(const char* name,
                const char* /*vol1Name*/, int /*vol1CopyNo*/,
                const char* /*vol2Name*/, int /*vol2CopyNo*/,
                const char* /*opSurfaceName*/) 
{
   Warning("SetBorderSurface",
           Form("Called for border surface %s. Not applicable in Geant3 - setting is ignored.", name));
}   
                
//______________________________________________________________________
void  TGeant3::SetSkinSurface(const char* name,
                const char* /*volName*/,
                const char* /*opSurfaceName*/) 
{
   Warning("SetSkinSurface",
           Form("Called for skin surface %s. Not applicable in Geant3 - setting is ignored.", name));
}   
                
//______________________________________________________________________
void  TGeant3::SetMaterialProperty(
                Int_t itmed, const char* /*propertyName*/, 
                Int_t /*np*/, Double_t* /*pp*/, Double_t* /*values*/) 
{
   Warning("SetMaterialProperty",
           Form("Called for material ID %5d. Not applicable in Geant3 - setting is ignored.", itmed));
}   
                
//______________________________________________________________________
void  TGeant3::SetMaterialProperty(
                Int_t itmed, const char* /*propertyName*/, 
                Double_t /*value*/) 
{
   Warning("SetMaterialProperty",
           Form("Called for material ID %5d. Not applicable in Geant3 - setting is ignored.", itmed));
}   
                
//______________________________________________________________________
void  TGeant3::SetMaterialProperty(
                const char* surfaceName, const char* /*propertyName*/, 
                Int_t /*np*/, Double_t* /*pp*/, Double_t* /*values*/) 
		{
   Warning("SetMaterialProperty",
           Form("Called for material surface  %s. Not applicable in Geant3 - setting is ignored.", surfaceName));
}   

//______________________________________________________________________
void  TGeant3::Gstpar(Int_t itmed, const char *param, Double_t parval)
{
  //
  //  To change the value of cut  or mechanism "CHPAR"
  //      to a new value PARVAL  for tracking medium ITMED
  //    The  data   structure  JTMED   contains  the   standard  tracking
  //  parameters (CUTS and flags to control the physics processes)  which
  //  are used  by default  for all  tracking media.   It is  possible to
  //  redefine individually  with GSTPAR  any of  these parameters  for a
  //  given tracking medium.
  //  ITMED     tracking medium number
  //  CHPAR     is a character string (variable name)
  //  PARVAL    must be given as a floating point.
  //

  Float_t fparval = parval;
  g3stpar(itmed,PASSCHARD(param), fparval PASSCHARL(param));
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gfkine(Int_t itra, Float_t *vert, Float_t *pvert, 
                      Int_t &ipart, Int_t &nvert)
{
  //           Storing/Retrieving Vertex and Track parameters
  //           ----------------------------------------------
  //
  //  Stores vertex parameters.
  //  VERT      array of (x,y,z) position of the vertex
  //  NTBEAM    beam track number origin of the vertex
  //            =0 if none exists
  //  NTTARG    target track number origin of the vertex
  //  UBUF      user array of NUBUF floating point numbers
  //  NUBUF
  //  NVTX      new vertex number (=0 in case of error).
  //  Prints vertex parameters.
  //  IVTX      for vertex IVTX.
  //            (For all vertices if IVTX=0)
  //  Stores long life track parameters.
  //  PLAB      components of momentum
  //  IPART     type of particle (see GSPART)
  //  NV        vertex number origin of track
  //  UBUF      array of NUBUF floating point user parameters
  //  NUBUF
  //  NT        track number (if=0 error).
  //  Retrieves long life track parameters.
  //  ITRA      track number for which parameters are requested
  //  VERT      vector origin of the track
  //  PVERT     4 momentum components at the track origin
  //  IPART     particle type (=0 if track ITRA does not exist)
  //  NVERT     vertex number origin of the track
  //  UBUF      user words stored in GSKINE.
  //  Prints initial track parameters.
  //  ITRA      for track ITRA
  //            (For all tracks if ITRA=0)
  //
  Float_t *ubuf=0;
  Int_t   nbuf;
  g3fkine(itra,vert,pvert,ipart,nvert,ubuf,nbuf);
}

//______________________________________________________________________
void  TGeant3::Gfvert(Int_t nvtx, Float_t *v, Int_t &ntbeam, Int_t &nttarg,
		      Float_t &tofg)
{
  //
  //       Retrieves the parameter of a vertex bank
  //       Vertex is generated from tracks NTBEAM NTTARG
  //       NVTX is the new vertex number
  //
  Float_t *ubuf=0;
  Int_t   nbuf;
  g3fvert(nvtx,v,ntbeam,nttarg,tofg,ubuf,nbuf);
}

//______________________________________________________________________
Int_t TGeant3::Gskine(Float_t *plab, Int_t ipart, Int_t nv, Float_t *buf,
		      Int_t nwbuf)
{
  //
  //       Store kinematics of track NT into data structure
  //       Track is coming from vertex NV
  //
  Int_t nt = 0;
  g3skine(plab, ipart, nv, buf, nwbuf, nt);
  return nt;
}

//______________________________________________________________________
Int_t TGeant3::Gsvert(Float_t *v, Int_t ntbeam, Int_t nttarg, Float_t *ubuf,
		      Int_t nwbuf)
{
  //
  //       Creates a new vertex bank
  //       Vertex is generated from tracks NTBEAM NTTARG
  //       NVTX is the new vertex number
  //
  Int_t nwtx = 0;
  g3svert(v, ntbeam, nttarg, ubuf, nwbuf, nwtx);
  return nwtx;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GPHYS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gphysi()
{
  //
  //       Initialize material constants for all the physics
  //       mechanisms used by GEANT
  //
  g3physi();
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GTRAK
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gdebug()
{
  //
  // Debug the current step
  //
  g3debug();
}

//______________________________________________________________________
void  TGeant3::Gekbin()
{
  //
  //       To find bin number in kinetic energy table
  //       stored in ELOW(NEKBIN)
  //
  g3ekbin();
}

//______________________________________________________________________
void  TGeant3::Gfinds()
{
  //
  //       Returns the set/volume parameters corresponding to
  //       the current space point in /GCTRAK/
  //       and fill common /GCSETS/
  //
  //       IHSET  user set identifier
  //       IHDET  user detector identifier
  //       ISET set number in JSET
  //       IDET   detector number in JS=LQ(JSET-ISET)
  //       IDTYPE detector type (1,2)
  //       NUMBV  detector volume numbers (array of length NVNAME)
  //       NVNAME number of volume levels
  //
  g3finds();
}

//______________________________________________________________________
void  TGeant3::Gsking(Int_t igk)
{
  //
  //   Stores in stack JSTAK either the IGKth track of /GCKING/,
  //    or the NGKINE tracks when IGK is 0.
  //
  g3sking(igk);
}

//______________________________________________________________________
void  TGeant3::Gskpho(Int_t igk)
{
  //
  //  Stores in stack JSTAK either the IGKth Cherenkov photon of
  //  /GCKIN2/, or the NPHOT tracks when IGK is 0.
  //
  g3skpho(igk);
}

//______________________________________________________________________
void  TGeant3::Gsstak(Int_t iflag)
{
  //
  //   Stores in auxiliary stack JSTAK the particle currently
  //    described in common /GCKINE/.
  //
  //   On request, creates also an entry in structure JKINE :
  //    IFLAG =
  //     0 : No entry in JKINE structure required (user)
  //     1 : New entry in JVERTX / JKINE structures required (user)
  //    <0 : New entry in JKINE structure at vertex -IFLAG (user)
  //     2 : Entry in JKINE structure exists already (from GTREVE)
  //
  g3sstak(iflag);
}

//______________________________________________________________________
void  TGeant3::Gsxyz()
{
  //
  //   Store space point VECT in banks JXYZ
  //
  g3sxyz();
}

//______________________________________________________________________
void  TGeant3::Gtrack()
{
  //
  //   Controls tracking of current particle
  //
  g3track();
}

//______________________________________________________________________
void  TGeant3::Gtreve()
{
  //
  //   Controls tracking of all particles belonging to the current event
  //
  g3treve();
}

//______________________________________________________________________
void  TGeant3::GtreveRoot()
{
  //
  //   Controls tracking of all particles belonging to the current event
  //
  gtreveroot();
}

//______________________________________________________________________
void  TGeant3::Grndm(Float_t *rvec, Int_t len) const
{
  //
  //  To set/retrieve the seed of the random number generator
  //
  TRandom* r=gMC->GetRandom();
  for(Int_t i=0; i<len; rvec[i++]=r->Rndm()) {};
}

//______________________________________________________________________
void  TGeant3::Grndmq(Int_t &is1, Int_t &is2, Int_t /*iseq*/,
		      const Text_t */*chopt*/)
{
  //
  //  To set/retrieve the seed of the random number generator
  //
  /*printf("Dummy grndmq called\n");*/
   is1 = gRandom->GetSeed();
   is2 = 0;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GDRAW
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gdxyz(Int_t /* it */)
{
  //
  // Draw the points stored with Gsxyz relative to track it
  //
}

//______________________________________________________________________
void  TGeant3::Gdcxyz()
{
  //
  // Draw the position of the current track
  //
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GGEOM
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void  TGeant3::Gdtom(Float_t *xd, Float_t *xm, Int_t iflag)
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
  //      IFLAG=2  convert direction cosines
  //
  g3dtom(xd, xm, iflag);
}

//______________________________________________________________________
void  TGeant3::Gdtom(Double_t *xd, Double_t *xm, Int_t iflag)
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
  //      IFLAG=2  convert direction cosines
  //

  Float_t* fxd = CreateFloatArray(xd, 3);
  Float_t* fxm = CreateFloatArray(xm, 3);

  Gdtom(fxd, fxm, iflag) ;

  for (Int_t i=0; i<3; i++) {
     xd[i] = fxd[i]; xm[i] = fxm[i];
  }

  delete [] fxd;
  delete [] fxm;
}

//______________________________________________________________________
void  TGeant3::Glmoth(const char* iudet, Int_t iunum, Int_t &nlev, 
                      Int_t *lvols, Int_t *lindx)
{
  //
  //   Loads the top part of the Volume tree in LVOLS (IVO's),
  //   LINDX (IN indices) for a given volume defined through
  //   its name IUDET and number IUNUM.
  //
  //   The routine stores only up to the last level where JVOLUM
  //   data structure is developed. If there is no development
  //   above the current level, it returns NLEV zero.
  Int_t *idum=0;
  g3lmoth(PASSCHARD(iudet), iunum, nlev, lvols, lindx, idum PASSCHARL(iudet));
}

//______________________________________________________________________
void  TGeant3::Gmedia(Float_t *x, Int_t &numed)
{
  //
  //   Finds in which volume/medium the point X is, and updates the
  //    common /GCVOLU/ and the structure JGPAR accordingly.
  //
  //   NUMED returns the tracking medium number, or 0 if point is
  //         outside the experimental setup.
  //

  static Int_t check = 0;
  g3media(x,numed,check);
}

//______________________________________________________________________
void  TGeant3::Gmtod(Float_t *xm, Float_t *xd, Int_t iflag)
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
  //           IFLAG=2  convert direction cosines
  //
  g3mtod(xm, xd, iflag);
}

//______________________________________________________________________
void  TGeant3::Gmtod(Double_t *xm, Double_t *xd, Int_t iflag)
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
  //           IFLAG=2  convert direction cosines
  //


  Float_t* fxm = CreateFloatArray(xm, 3);
  Float_t* fxd = CreateFloatArray(xd, 3);

  Gmtod(fxm, fxd, iflag) ;

  for (Int_t i=0; i<3; i++) {
     xm[i] = fxm[i]; xd[i] = fxd[i];
  }

  delete [] fxm;
  delete [] fxd;
}

//______________________________________________________________________
void  TGeant3::Gsdvn(const char *name, const char *mother, Int_t ndiv,
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
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  g3sdvn(PASSCHARD(vname), PASSCHARD(vmother), ndiv, iaxis PASSCHARL(vname)
	PASSCHARL(vmother));
}

//______________________________________________________________________
void  TGeant3::Gsdvn2(const char *name, const char *mother, Int_t ndiv,
		      Int_t iaxis, Double_t c0i, Int_t numed)
{
  //
  // Create a new volume by dividing an existing one
  //
  // Divides mother into ndiv divisions called name
  // along axis iaxis starting at coordinate value c0.
  // the new volume created will be medium number numed.
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  Float_t fc0i = c0i;
  g3sdvn2(PASSCHARD(vname), PASSCHARD(vmother), ndiv, iaxis, fc0i, numed
	 PASSCHARL(vname) PASSCHARL(vmother));
}

//______________________________________________________________________
void  TGeant3::Gsdvs(const char *name, const char *mother, Float_t step,
		     Int_t iaxis, Int_t numed)
{
  //
  // Create a new volume by dividing an existing one
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  g3sdvs(PASSCHARD(vname), PASSCHARD(vmother), step, iaxis, numed
	PASSCHARL(vname) PASSCHARL(vmother));
}

//______________________________________________________________________
void  TGeant3::Gsdvs2(const char *name, const char *mother, Float_t step,
		      Int_t iaxis, Float_t c0, Int_t numed)
{
  //
  // Create a new volume by dividing an existing one
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  g3sdvs2(PASSCHARD(vname), PASSCHARD(vmother), step, iaxis, c0, numed
	 PASSCHARL(vname) PASSCHARL(vmother));
}

//______________________________________________________________________
void  TGeant3::Gsdvt(const char *name, const char *mother, Double_t step,
		     Int_t iaxis, Int_t numed, Int_t ndvmx)
{
  //
  // Create a new volume by dividing an existing one
  //
  //       Divides MOTHER into divisions called NAME along
  //       axis IAXIS in steps of STEP. If not exactly divisible
  //       will make as many as possible and will center them
  //       with respect to the mother. Divisions will have medium
  //       number NUMED. If NUMED is 0, NUMED of MOTHER is taken.
  //       NDVMX is the expected maximum number of divisions
  //          (If 0, no protection tests are performed)
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  Float_t fstep = step;
  g3sdvt(PASSCHARD(vname), PASSCHARD(vmother), fstep, iaxis, numed, ndvmx
	PASSCHARL(vname) PASSCHARL(vmother));
}

//______________________________________________________________________
void  TGeant3::Gsdvt2(const char *name, const char *mother, Double_t step,
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
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  Float_t fstep = step;
  Float_t fc0 = c0;
  g3sdvt2(PASSCHARD(vname), PASSCHARD(vmother), fstep, iaxis, fc0,
	 numed, ndvmx PASSCHARL(vname) PASSCHARL(vmother));
}

//______________________________________________________________________
void  TGeant3::Gsord(const char *name, Int_t iax)
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

  char vname[5];
  Vname(name,vname);
  g3sord(PASSCHARD(vname), iax PASSCHARL(vname));
}

//______________________________________________________________________
void  TGeant3::Gspos(const char *name, Int_t nr, const char *mother, 
                     Double_t x, Double_t y, Double_t z, Int_t irot, 
                     const char *konly)
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

  TString only = konly;
  only.ToLower();
  Bool_t isOnly = kFALSE;
  if (only.Contains("only")) isOnly = kTRUE;
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  Float_t fx = x;
  Float_t fy = y;
  Float_t fz = z;
  g3spos(PASSCHARD(vname), nr, PASSCHARD(vmother), fx, fy, fz, irot,
	PASSCHARD(konly) PASSCHARL(vname) PASSCHARL(vmother)
	PASSCHARL(konly));
}

//______________________________________________________________________
void  TGeant3::G3Gsposp(const char *name, Int_t nr, const char *mother,
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Float_t *upar, Int_t np )
{
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //
  TString only = konly;
  only.ToLower();
  Bool_t isOnly = kFALSE;
  if (only.Contains("only")) isOnly = kTRUE;
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  Float_t fx = x;
  Float_t fy = y;
  Float_t fz = z;
  g3sposp(PASSCHARD(vname), nr, PASSCHARD(vmother), fx, fy, fz, irot,
	 PASSCHARD(konly), upar, np PASSCHARL(vname) PASSCHARL(vmother)
	 PASSCHARL(konly));
}

//______________________________________________________________________
void  TGeant3::Gsposp(const char *name, Int_t nr, const char *mother,
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Float_t *upar, Int_t np )
{
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //

  G3Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np);
}

//______________________________________________________________________
void  TGeant3::Gsposp(const char *name, Int_t nr, const char *mother,
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Double_t *upar, Int_t np )
{
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //

  Float_t* fupar = CreateFloatArray(upar, np);
  G3Gsposp(name, nr, mother, x, y, z, irot, konly, fupar, np);
  delete [] fupar;
}

//______________________________________________________________________
void  TGeant3::Gsrotm(Int_t nmat, Float_t theta1, Float_t phi1, Float_t theta2,
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

  g3srotm(nmat, theta1, phi1, theta2, phi2, theta3, phi3);
}

//______________________________________________________________________
void  TGeant3::Gprotm(Int_t nmat)
{
  //
  //    To print rotation matrices structure JROTM
  //     nmat     Rotation matrix number
  //
  g3protm(nmat);
 }

//______________________________________________________________________
Int_t TGeant3::G3Gsvolu(const char *name, const char *shape, Int_t nmed,
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
  char vname[5];
  Vname(name,vname);
  char vshape[5];
  Vname(shape,vshape);

  g3svolu(PASSCHARD(vname), PASSCHARD(vshape), nmed, upar, npar, ivolu
	 PASSCHARL(vname) PASSCHARL(vshape));

  return ivolu;
}
//______________________________________________________________________
Int_t TGeant3::Gsvolu(const char *name, const char *shape, Int_t nmed,
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
  ivolu = G3Gsvolu(name, shape, nmed, upar, npar);
  return ivolu;

}

//______________________________________________________________________
Int_t TGeant3::Gsvolu(const char *name, const char *shape, Int_t nmed,
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
  Float_t* fupar = CreateFloatArray(upar, npar);
  ivolu = G3Gsvolu(name, shape, nmed, fupar, npar);
  delete [] fupar;
  return ivolu;
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//           T H E    D R A W I N G   P A C K A G E
//           ======================================
//  Drawing functions. These functions allow the visualization in several 
//  ways of the volumes defined in the geometrical data structure. It is 
//  possible to draw the logical tree of volumes belonging to the detector 
//  (DTREE), to show their geometrical specification (DSPEC,DFSPC), to 
//  draw them and their cut views (DRAW, DCUT). Moreover, it is possible 
//  to execute these commands when the hidden line removal option is 
//  activated; in this case, the volumes can be also either translated 
//  in the space (SHIFT), or clipped by boolean operation (CVOL). In 
//  addition, it is possible to fill the surfaces of the volumes
//  with solid colors when the shading option (SHAD) is activated.
//  Several tools (ZOOM, LENS) have been developed to zoom detailed parts
//  of the detectors or to scan physical events as well.
//  Finally, the command MOVE will allow the rotation, translation and 
//  zooming on real time parts of the detectors or tracks and hits of a 
//  simulated event. Ray-tracing commands. In case the command (DOPT RAYT 
//  ON) is executed, the drawing is performed by the Geant ray-tracing;
//  automatically, the color is assigned according to the tracking medium 
//  of each volume and the volumes with a density lower/equal than the 
//  air are considered transparent; if the option (USER) is set (ON) 
//  (again via the command (DOPT)), the user can set color and visibility 
//  for the desired volumes via the command (SATT), as usual, relatively 
//  to the attributes (COLO) and (SEEN). The resolution can be set via 
//  the command (SATT * FILL VALUE), where (VALUE) is the ratio between 
//  the number of pixels drawn and 20 (user coordinates). Parallel view 
//  and perspective view are possible (DOPT PROJ PARA/PERS); in the first 
//  case, we assume that the first mother volume of the tree is a box with
//  dimensions 10000 X 10000 X 10000 cm and the view point (infinitely far) 
//  is 5000 cm far from the origin along the Z axis of the user coordinates; 
//  in the second case, the distance between the observer and the origin 
//  of the world reference system is set in cm by the command (PERSP NAME 
//  VALUE); grand-angle or telescopic effects can be achieved changing the 
//  scale factors in the command (DRAW). When the final picture does not 
//  occupy the full window, mapping the space before tracing can speed up 
//  the drawing, but can also produce less precise results; values from 1 
//  to 4 are allowed in the command (DOPT MAPP VALUE), the mapping being 
//  more precise for increasing (VALUE); for (VALUE = 0) no mapping is 
//  performed (therefore max precision and lowest speed). The command 
//  (VALCUT) allows the cutting of the detector by three planes orthogonal 
//  to the x,y,z axis. The attribute (LSTY) can be set by the command
//  SATT for any desired volume and can assume values from 0 to 7; it 
//  determines the different light processing to be performed for different 
//  materials:
//  0 = dark-matt, 1 = bright-matt, 2 = plastic, 3 = ceramic, 4 = rough-metals,
//  5 = shiny-metals, 6 = glass, 7 = mirror. The detector is assumed to 
//  be in the dark, the ambient light luminosity is 0.2 for each basic 
//  hue (the saturation is 0.9) and the observer is assumed to have a 
//  light source (therefore he will produce parallel light in the case 
//  of parallel view and point-like-source light in the case of perspective 
//  view).
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________
void TGeant3::Gsatt(const char *name, const char *att, Int_t val)
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
  //            represent the line width of the scan lines filling the surfaces
  //            (whereas the FILL value represent their number). Therefore
  //            tuning this parameter will help to obtain the desired
  //            quality/performance ratio.
  //
  //     COLO   color code -166,...,1,2,..166 (default=1)
  //            n=1=black
  //            n=2=red;    n=17+m, m=0,25, increasing luminosity according to 'm';
  //            n=3=green;  n=67+m, m=0,25, increasing luminosity according to 'm';
  //            n=4=blue;   n=117+m, m=0,25, increasing luminosity according to 'm';
  //            n=5=yellow; n=42+m, m=0,25, increasing luminosity according to 'm';
  //            n=6=violet; n=142+m, m=0,25, increasing luminosity according to 'm';
  //            n=7=light-blue; n=92+m, m=0,25, increasing luminosity according to 'm';
  //            color=n*10+m, m=1,2,...9, will produce the same color
  //            as 'n', but with increasing luminosity according to 'm';
  //            COLO<0 will act as if abs(COLO) was set for the volume
  //            and for all the levels below it.
  //            When for a volume the attribute FILL is > 1 (and the
  //            option SHAD is on), the ABS of its color code must be < 8
  //            because an automatic shading of its faces will be
  //            performed.
  //
  //     FILL  (1992) fill area  -7,...,0,1,...7 (default=0)
  //            when option SHAD is "on" the FILL attribute of any
  //            volume can be set different from 0 (normal drawing);
  //            if it is set to 1, the faces of such volume will be filled
  //            with solid colors; if ABS(FILL) is > 1, then a light
  //            source is placed along the observer line, and the faces of
  //            such volumes will be painted by colors whose luminosity
  //            will depend on the amount of light reflected;
  //            if ABS(FILL) = 1, then it is possible to use all the 166
  //            colors of the color table, because the automatic shading
  //            is not performed;
  //            for increasing values of FILL the drawing will be performed
  //            with higher and higher resolution improving the quality (the
  //            number of scan lines used to fill the faces increases with 
  //            FILL); it is possible to set different values of FILL
  //            for different volumes, in order to optimize at the same time
  //            the performance and the quality of the picture;
  //            FILL<0 will act as if abs(FILL) was set for the volume
  //            and for all the levels below it.
  //            This kind of drawing can be saved in 'picture files'
  //            or in view banks.
  //            0=drawing without fill area
  //            1=faces filled with solid colors and resolution = 6
  //            2=lowest resolution (very fast)
  //            3=default resolution
  //            4=.................
  //            5=.................
  //            6=.................
  //            7=max resolution
  //            Finally, if a colored background is desired, the FILL
  //            attribute for the first volume of the tree must be set
  //            equal to -abs(colo), colo being >0 and <166.
  //
  //     SET   set number associated to volume name
  //     DET   detector number associated to volume name
  //     DTYP  detector type (1,2)
  //

  char vname[5];
  Vname(name,vname);
  char vatt[5];
  Vname(att,vatt);
  g3satt(PASSCHARD(vname), PASSCHARD(vatt), val PASSCHARL(vname)
	PASSCHARL(vatt));
}

//______________________________________________________________________
void TGeant3::Gfpara(const char *name, Int_t number, Int_t intext, Int_t& npar,
			 Int_t& natt, Float_t* par, Float_t* att)
{
  //
  // Find the parameters of a volume
  //
  g3fpara(PASSCHARD(name), number, intext, npar, natt, par, att
	 PASSCHARL(name));
}

//______________________________________________________________________
void TGeant3::Gckpar(Int_t ish, Int_t npar, Float_t* par)
{
  //
  // Check the parameters of a shape
  //
  gckpar(ish,npar,par);
}

//______________________________________________________________________
void TGeant3::Gckmat(Int_t itmed, char* natmed)
{
  //
  // Check the parameters of a tracking medium
  //
  g3ckmat(itmed, PASSCHARD(natmed) PASSCHARL(natmed));
}

//______________________________________________________________________
Int_t TGeant3::Glvolu(Int_t nlev, Int_t *lnam,Int_t *lnum)
{
  //
  //  nlev   number of levels deep into the volume tree
  //         size of the arrays lnam and lnum
  //  lnam   an integer array who's 4 bytes contain the ASCII code for the
  //         volume names
  //  lnum   an integer array containing the copy numbers for that specific
  //         volume
  //
  //  This routine fills the volume parameters in common /gcvolu/ for a
  //  physical tree, specified by the list lnam and lnum of volume names
  //  and numbers, and for all its ascendants up to level 1. This routine
  //  is optimized and does not re-compute the part of the history already
  //  available in GCVOLU. This means that if it is used in user programs
  //  outside the usual framework of the tracking, the user has to initialize
  //  to zero NLEVEL in the common GCVOLU. It return 0 if there were no
  //  problems in make the call.
  //
  Int_t ier;
  g3lvolu(nlev, lnam, lnum, ier);
  return ier;
}

//______________________________________________________________________
void TGeant3::Gdelete(Int_t /* iview */)
{
  //
  //  IVIEW  View number
  //
  //  It deletes a view bank from memory.
  //
}

//______________________________________________________________________
void TGeant3::Gdopen(Int_t /* iview */)
{
  //
  //  IVIEW  View number
  //
  //  When a drawing is very complex and requires a long time to be
  //  executed, it can be useful to store it in a view bank: after a
  //  call to DOPEN and the execution of the drawing (nothing will
  //  appear on the screen), and after a necessary call to DCLOSE,
  //  the contents of the bank can be displayed in a very fast way
  //  through a call to DSHOW; therefore, the detector can be easily
  //  zoomed many times in different ways. Please note that the pictures
  //  with solid colors can now be stored in a view bank or in 'PICTURE FILES'
  //
}

//______________________________________________________________________
void TGeant3::Gdclose()
{
  //
  //  It closes the currently open view bank; it must be called after the
  //  end of the drawing to be stored.
  //
}

//______________________________________________________________________
void TGeant3::Gdshow(Int_t /* iview */)
{
  //
  //  IVIEW  View number
  //
  //  It shows on the screen the contents of a view bank. It
  //  can be called after a view bank has been closed.
  //
}

//______________________________________________________________________
void TGeant3::Gdopt(const char *name,const char *value)
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
  //     USER    ON       User graphics options in the ray tracing.
  //             OFF (D)  Automatic graphics options.
  //

  char vname[5];
  Vname(name,vname);
  char vvalue[5];
  Vname(value,vvalue);
  //g3dopt(PASSCHARD(vname), PASSCHARD(vvalue) PASSCHARL(vname)
  //	PASSCHARL(vvalue));
}

//______________________________________________________________________
void TGeant3::Gdraw(const char* /*name*/,Double_t /*theta*/, Double_t /*phi*/, 
                    Double_t /*psi*/, Double_t /*u0*/, Double_t /*v0*/, Double_t /*ul*/,
                    Double_t /*vl*/)
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
  //  size of its ZEBRA store). Finally, at the end of each drawing (with 
  //  HIDE on), the program will print messages about the memory used and
  //  statistics on the volumes' visibility.
  //  The following commands will produce the drawing of a green
  //  volume, specified by NAME, without using the hidden line removal
  //  technique, using the hidden line removal technique,
  //  with different line width and color (red), with
  //  solid color, with shading of surfaces, and without edges.
  //  Finally, some examples are given for the ray-tracing. (A possible
  //  string for the NAME of the volume can be found using the command DTREE).
  //
}

//______________________________________________________________________
void TGeant3::Gdrawc(const char* /*name*/,Int_t /*axis*/, Float_t /*cut*/, Float_t /*u0*/,
		     Float_t /*v0*/, Float_t /*ul*/, Float_t /*vl*/)
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
}

//______________________________________________________________________
void TGeant3::Gdrawx(const char* /*name*/, Float_t /*cutthe*/, Float_t /*cutphi*/,
		     Float_t /*cutval*/, Float_t /*theta*/, Float_t /*phi*/, Float_t /*u0*/,
		     Float_t /*v0*/,Float_t /*ul*/, Float_t /*vl*/)
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
}

//______________________________________________________________________
void TGeant3::Gdhead(Int_t /*isel*/, const char* /*name*/, Double_t /*chrsiz*/)
{
  //
  //  Parameters
  //  +
  //  ISEL   Option flag  D=111110
  //  NAME   Title
  //  CHRSIZ Character size (cm) of title NAME D=0.6
  //
  //  ISEL =
  //   0      to have only the header lines
  //   xxxxx1 to add the text name centered on top of header
  //   xxxx1x to add global detector name (first volume) on left
  //   xxx1xx to add date on right
  //   xx1xxx to select thick characters for text on top of header
  //   x1xxxx to add the text 'EVENT NR x' on top of header
  //   1xxxxx to add the text 'RUN NR x' on top of header
  //  NOTE that ISEL=x1xxx1 or ISEL=1xxxx1 are illegal choices,
  //  i.e. they generate overwritten text.
  //
}

//______________________________________________________________________
void TGeant3::Gdman(Double_t /*u*/, Double_t /*v*/, const char* /*type*/)
{
  //
  //  Draw a 2D-man at position (U0,V0)
  //  Parameters
  //  U      U-coord. (horizontal) of the center of man' R
  //  V      V-coord. (vertical) of the center of man' R
  //  TYPE   D='MAN' possible values: 'MAN,WM1,WM2,WM3'
  //
  //   CALL GDMAN(u,v),CALL GDWMN1(u,v),CALL GDWMN2(u,v),CALL GDWMN2(u,v)
  //  It superimposes the picture of a man or of a woman, chosen among
  //  three different ones, with the same scale factors as the detector
  //  in the current drawing.
  //
}

//______________________________________________________________________
void TGeant3::Gdspec(const char* /*name*/)
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
}

//______________________________________________________________________
void TGeant3::DrawOneSpec(const char* /*name*/)
{
  //
  //  Function called when one double-clicks on a volume name
  //  in a TPavelabel drawn by Gdtree.
  //
}

//______________________________________________________________________
void TGeant3::Gdtree(const char* /*name*/, Int_t /*levmax*/, Int_t /*isel*/)
{
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree,
  //  Each volume in the tree is represented by a TPaveTree object.
  //  Double-clicking on a TPaveTree draws the specs of the corresponding 
  //  volume.
  //  Use TPaveTree pop-up menu to select:
  //    - drawing specs
  //    - drawing tree
  //    - drawing tree of parent
  //
}

//______________________________________________________________________
void TGeant3::GdtreeParent(const char* /*name*/, Int_t /*levmax*/, Int_t /*isel*/)
{
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree of the parent of name.
  //
}

//______________________________________________________________________
void TGeant3::SetABAN(Int_t par)
{
  //
  // par = 1 particles will be stopped according to their residual
  //         range if they are not in a sensitive material and are
  //         far enough from the boundary
  //       0 particles are transported normally
  //
  fGcphys->dphys1 = par;
  SetBit(kABAN);
}


//______________________________________________________________________
void TGeant3::SetANNI(Int_t par)
{
  //
  //   To control positron annihilation.
  //    par =0 no annihilation
  //        =1 annihilation. Decays processed.
  //        =2 annihilation. No decay products stored.
  //
  fGcphys->ianni = par;
}


//______________________________________________________________________
void TGeant3::SetAUTO(Int_t par)
{
  //
  //  To control automatic calculation of tracking medium parameters:
  //   par =0 no automatic calculation;
  //       =1 automatic calculation.
  //
  fGctrak->igauto = par;
  SetBit(kAUTO);
}


//______________________________________________________________________
void TGeant3::SetBOMB(Float_t /*boom*/)
{
  //
  //  BOOM  : Exploding factor for volumes position
  //
  //  To 'explode' the detector. If BOOM is positive (values smaller
  //  than 1. are suggested, but any value is possible)
  //  all the volumes are shifted by a distance
  //  proportional to BOOM along the direction between their center
  //  and the origin of the MARS; the volumes which are symmetric
  //  with respect to this origin are simply not shown.
  //  BOOM equal to 0 resets the normal mode.
  //  A negative (greater than -1.) value of
  //  BOOM will cause an 'implosion'; for even lower values of BOOM
  //  the volumes' positions will be reflected respect to the origin.
  //  This command can be useful to improve the 3D effect for very
  //  complex detectors. The following commands will make explode the
  //  detector:
  //
}

//______________________________________________________________________
void TGeant3::SetBREM(Int_t par)
{
  //
  //  To control bremsstrahlung.
  //   par =0 no bremsstrahlung
  //       =1 bremsstrahlung. Photon processed.
  //       =2 bremsstrahlung. No photon stored.
  //
  fGcphys->ibrem = par;
}


//______________________________________________________________________
void TGeant3::SetCKOV(Int_t par)
{
  //
  //  To control Cerenkov production
  //   par =0 no Cerenkov;
  //       =1 Cerenkov;
  //       =2 Cerenkov with primary stopped at each step.
  //
  fGctlit->itckov = par;
}


//______________________________________________________________________
void  TGeant3::SetClipBox(const char* /*name*/, Double_t /*xmin*/, Double_t /*xmax*/,
			  Double_t /*ymin*/, Double_t /*ymax*/, Double_t /*zmin*/, Double_t /*zmax*/)
{
  //
  //  The hidden line removal technique is necessary to visualize properly
  //  very complex detectors. At the same time, it can be useful to visualize
  //  the inner elements of a detector in detail. This function allows
  //  subtractions (via boolean operation) of BOX shape from any part of
  //  the detector, therefore showing its inner contents.
  //  If "*" is given as the name of the
  //  volume to be clipped, all volumes are clipped by the given box.
  //  A volume can be clipped at most twice.
  //  if a volume is explicitly clipped twice,
  //  the "*" will not act on it anymore. Giving "." as the name
  //  of the volume to be clipped will reset the clipping.
  //  Parameters
  //  NAME   Name of volume to be clipped
  //  +
  //  XMIN   Lower limit of the Shape X coordinate
  //  XMAX   Upper limit of the Shape X coordinate
  //  YMIN   Lower limit of the Shape Y coordinate
  //  YMAX   Upper limit of the Shape Y coordinate
  //  ZMIN   Lower limit of the Shape Z coordinate
  //  ZMAX   Upper limit of the Shape Z coordinate
  //
  //  This function performs a boolean subtraction between the volume
  //  NAME and a box placed in the MARS according the values of the given
  //  coordinates.

}

//______________________________________________________________________
void TGeant3::SetCOMP(Int_t par)
{
  //
  //  To control Compton scattering
  //   par =0 no Compton
  //       =1 Compton. Electron processed.
  //       =2 Compton. No electron stored.
  //
  //
  fGcphys->icomp = par;
}

//modified by Andrea Fontana and Alberto Rotondi - march 2007
//added array of 5 user definable cuts (like in old Geant)
//______________________________________________________________________
void TGeant3::SetCUTS(Float_t cutgam,Float_t cutele,Float_t cutneu,
		      Float_t cuthad,Float_t cutmuo ,Float_t bcute ,
		      Float_t bcutm ,Float_t dcute ,Float_t dcutm ,
		      Float_t ppcutm, Float_t tofmax, Float_t *gcuts)
{
  //
  //  CUTGAM   Cut for gammas              D=0.001
  //  CUTELE   Cut for electrons           D=0.001
  //  CUTHAD   Cut for charged hadrons     D=0.01
  //  CUTNEU   Cut for neutral hadrons     D=0.01
  //  CUTMUO   Cut for muons               D=0.01
  //  BCUTE    Cut for electron brems.     D=-1.
  //  BCUTM    Cut for muon brems.         D=-1.
  //  DCUTE    Cut for electron delta-rays D=-1.
  //  DCUTM    Cut for muon delta-rays     D=-1.
  //  PPCUTM   Cut for e+e- pairs by muons D=0.01
  //  TOFMAX   Time of flight cut          D=1.E+10
  //
  //   If the default values (-1.) for       BCUTE ,BCUTM ,DCUTE ,DCUTM
  //   are not modified, they will be set to CUTGAM,CUTGAM,CUTELE,CUTELE
  //   respectively.
  //  If one of the parameters from CUTGAM to PPCUTM included
  //  is modified, cross-sections and energy loss tables must be
  //  recomputed via the function Gphysi.
  //
  fGccuts->cutgam = cutgam;
  fGccuts->cutele = cutele;
  fGccuts->cutneu = cutneu;
  fGccuts->cuthad = cuthad;
  fGccuts->cutmuo = cutmuo;
  fGccuts->bcute  = bcute;
  fGccuts->bcutm  = bcutm;
  fGccuts->dcute  = dcute;
  fGccuts->dcutm  = dcutm;
  fGccuts->ppcutm = ppcutm;
  fGccuts->tofmax = tofmax;
  fGccuts->gcuts[0] = gcuts[0];
  fGccuts->gcuts[1] = gcuts[1];
  fGccuts->gcuts[2] = gcuts[2];
  fGccuts->gcuts[3] = gcuts[3];
  fGccuts->gcuts[4] = gcuts[4];
}

//added by Andrea Fontana and Alberto Rotondi - april 2007
//______________________________________________________________________
void TGeant3::SetECut(Float_t gcalpha)
{
  fGcmore->gcalpha = gcalpha;
}

void TGeant3::SetClose(Int_t iclose,Float_t *pf,Float_t dstrt,
                       Float_t *w1, Float_t *w2,
		       Float_t *p1,Float_t *p2,Float_t *p3,Float_t *clen)
{
  fGcmore->iclose = iclose;
  fGcmore->pfinal[0] = pf[0];
  fGcmore->pfinal[1] = pf[1];
  fGcmore->pfinal[2] = pf[2];
  fGcmore->dstrt = dstrt;
  fGcmore->wire1[0] = w1[0];
  fGcmore->wire1[1] = w1[1];
  fGcmore->wire1[2] = w1[2];
  fGcmore->wire2[0] = w2[0];
  fGcmore->wire2[1] = w2[1];
  fGcmore->wire2[2] = w2[2];
  fGcmore->p1[0] = p1[0];
  fGcmore->p1[1] = p1[1];
  fGcmore->p1[2] = p1[2];
  fGcmore->p2[0] = p2[0];
  fGcmore->p2[1] = p2[1];
  fGcmore->p2[2] = p2[2];
  fGcmore->p3[0] = p3[0];
  fGcmore->p3[1] = p3[1];
  fGcmore->p3[2] = p3[2];
  fGcmore->cleng[0] = clen[0];
  fGcmore->cleng[1] = clen[1];
  fGcmore->cleng[2] = clen[2];
}

void TGeant3::GetClose(Float_t *p1,Float_t *p2,Float_t *p3,Float_t *len)
{
  p1[0] = fGcmore->p1[0];
  p1[1] = fGcmore->p1[1];
  p1[2] = fGcmore->p1[2];
  p2[0] = fGcmore->p2[0];
  p2[1] = fGcmore->p2[1];
  p2[2] = fGcmore->p2[2];
  p3[0] = fGcmore->p3[0];
  p3[1] = fGcmore->p3[1];
  p3[2] = fGcmore->p3[2];
  len[0] = fGcmore->cleng[0];
  len[1] = fGcmore->cleng[1];
  len[2] = fGcmore->cleng[2];
}

//______________________________________________________________________
void TGeant3::SetDCAY(Int_t par)
{
  //
  //  To control Decay mechanism.
  //   par =0 no decays.
  //       =1 Decays. secondaries processed.
  //       =2 Decays. No secondaries stored.
  //
  fGcphys->idcay = par;
}


//______________________________________________________________________
void TGeant3::SetDEBU(Int_t emin, Int_t emax, Int_t emod)
{
  //
  // Set the debug flag and frequency
  // Selected debug output will be printed from
  // event emin to even emax each emod event
  //
  fGcflag->idemin = emin;
  fGcflag->idemax = emax;
  fGcflag->itest  = emod;
  SetBit(kDEBU);
}


//______________________________________________________________________
void TGeant3::SetDRAY(Int_t par)
{
  //
  //  To control delta rays mechanism.
  //   par =0 no delta rays.
  //       =1 Delta rays. secondaries processed.
  //       =2 Delta rays. No secondaries stored.
  //
  fGcphys->idray = par;
}

//______________________________________________________________________
void TGeant3::SetERAN(Float_t ekmin, Float_t ekmax, Int_t nekbin)
{
  //
  //  To control cross section tabulations
  //   ekmin = minimum kinetic energy in GeV
  //   ekmax = maximum kinetic energy in GeV
  //   nekbin = number of logarithmic bins (<200)
  //
  fGcmulo->ekmin = ekmin;
  fGcmulo->ekmax = ekmax;
  fGcmulo->nekbin = nekbin;
  SetBit(kERAN);
}

//______________________________________________________________________
void TGeant3::SetHADR(Int_t par)
{
  //
  //  To control hadronic interactions.
  //   par =0 no hadronic interactions.
  //       =1 Hadronic interactions. secondaries processed.
  //       =2 Hadronic interactions. No secondaries stored.
  //
  fGcphys->ihadr = par;
}

//______________________________________________________________________
void TGeant3::SetKINE(Int_t kine, Float_t xk1, Float_t xk2, Float_t xk3,
		      Float_t xk4, Float_t xk5, Float_t xk6, Float_t xk7,
		      Float_t xk8, Float_t xk9, Float_t xk10)
{
  //
  // Set the variables in /GCFLAG/ IKINE, PKINE(10)
  // Their meaning is user defined
  //
  fGckine->ikine    = kine;
  fGckine->pkine[0] = xk1;
  fGckine->pkine[1] = xk2;
  fGckine->pkine[2] = xk3;
  fGckine->pkine[3] = xk4;
  fGckine->pkine[4] = xk5;
  fGckine->pkine[5] = xk6;
  fGckine->pkine[6] = xk7;
  fGckine->pkine[7] = xk8;
  fGckine->pkine[8] = xk9;
  fGckine->pkine[9] = xk10;
}

//______________________________________________________________________
void TGeant3::SetLOSS(Int_t par)
{
  //
  //  To control energy loss.
  //   par =0 no energy loss;
  //       =1 restricted energy loss fluctuations;
  //       =2 complete energy loss fluctuations;
  //       =3 same as 1;
  //       =4 no energy loss fluctuations.
  //  If the value ILOSS is changed, then cross-sections and energy loss
  //  tables must be recomputed via the command 'PHYSI'.
  //
  fGcphys->iloss = par;
}


//______________________________________________________________________
void TGeant3::SetMULS(Int_t par)
{
  //
  //  To control multiple scattering.
  //   par =0 no multiple scattering.
  //       =1 Moliere or Coulomb scattering.
  //       =2 Moliere or Coulomb scattering.
  //       =3 Gaussian scattering.
  //
  fGcphys->imuls = par;
}


//______________________________________________________________________
void TGeant3::SetMUNU(Int_t par)
{
  //
  //  To control muon nuclear interactions.
  //   par =0 no muon-nuclear interactions.
  //       =1 Nuclear interactions. Secondaries processed.
  //       =2 Nuclear interactions. Secondaries not processed.
  //
  fGcphys->imunu = par;
}

//______________________________________________________________________
void TGeant3::SetOPTI(Int_t par)
{
  //
  //  This flag controls the tracking optimization performed via the
  //  GSORD routine:
  //      1 no optimization at all; GSORD calls disabled;
  //      0 no optimization; only user calls to GSORD kept;
  //      1 all non-GSORDered volumes are ordered along the best axis;
  //      2 all volumes are ordered along the best axis.
  //
  fGcopti->ioptim = par;
  SetBit(kOPTI);
}

//______________________________________________________________________
void TGeant3::SetPAIR(Int_t par)
{
  //
  //  To control pair production mechanism.
  //   par =0 no pair production.
  //       =1 Pair production. secondaries processed.
  //       =2 Pair production. No secondaries stored.
  //
  fGcphys->ipair = par;
}


//______________________________________________________________________
void TGeant3::SetPFIS(Int_t par)
{
  //
  //  To control photo fission mechanism.
  //   par =0 no photo fission.
  //       =1 Photo fission. secondaries processed.
  //       =2 Photo fission. No secondaries stored.
  //
  fGcphys->ipfis = par;
}

//______________________________________________________________________
void TGeant3::SetPHOT(Int_t par)
{
  //
  //  To control Photo effect.
  //   par =0 no photo electric effect.
  //       =1 Photo effect. Electron processed.
  //       =2 Photo effect. No electron stored.
  //
  fGcphys->iphot = par;
}

//______________________________________________________________________
void TGeant3::SetRAYL(Int_t par)
{
  //
  //  To control Rayleigh scattering.
  //   par =0 no Rayleigh scattering.
  //       =1 Rayleigh.
  //
  fGcphys->irayl = par;
}

//______________________________________________________________________
void TGeant3::SetSTRA(Int_t par)
{
  //
  //  To control energy loss fluctuations
  //  with the Photo-Absorption Ionization model.
  //   par =0 no Straggling.
  //       =1 Straggling yes => no Delta rays.
  //
  fGcphlt->istra = par;
}

//______________________________________________________________________
void TGeant3::SetSWIT(Int_t sw, Int_t val)
{
  //
  //  sw    Switch number
  //  val   New switch value
  //
  //  Change one element of array ISWIT(10) in /GCFLAG/
  //
  if (sw <= 0 || sw > 10) return;
  fGcflag->iswit[sw-1] = val;
  SetBit(kSWIT);
}


//______________________________________________________________________
void TGeant3::SetTRIG(Int_t nevents)
{
  //
  // Set number of events to be run
  //
  fGcflag->nevent = nevents;
  SetBit(kTRIG);
}

//______________________________________________________________________
void TGeant3::SetUserDecay(Int_t pdg)
{
  //
  // Force the decays of particles to be done with Pythia
  // and not with the Geant routines.
  // just kill pointers doing mzdrop
  //
  Int_t ipart = IdFromPDG(pdg);
  if(ipart<0) {
    printf("Particle %d not in geant\n",pdg);
    return;
  }
  Int_t jpart=fGclink->jpart;
  Int_t jpa=fZlq[jpart-ipart];
  //
  if(jpart && jpa) {
    Int_t jpa1=fZlq[jpa-1];
    if(jpa1)
      mzdrop(fGcbank->ixcons,jpa1,PASSCHARD(" ") PASSCHARL(" "));
    Int_t jpa2=fZlq[jpa-2];
    if(jpa2)
      mzdrop(fGcbank->ixcons,jpa2,PASSCHARD(" ") PASSCHARL(" "));
  }
}
//______________________________________________________________________
Bool_t TGeant3::SetDecayMode(Int_t pdg, Float_t bratio[6], Int_t mode[6][3])
{
  //
  // Set user decay modes by calling Gsdk
  //
   if ( pdg == 0 ) {
     printf("Cannot define decay mode for particle with PDG=0");
     return false;
   }

   if ( IdFromPDG(pdg) < 0 ) {
     printf("Particle %d not in geant\n",pdg);
     return false;
   }

   SetUserDecay(pdg);

   Int_t g3mode[6];
   Int_t id1,id2,id3;
   for (Int_t k1=0; k1<6; k1++) g3mode[k1]=0;
   for (Int_t k=0; k<6; k++) {

      if(mode[k][0]!=0) {
        id1= IdFromPDG(mode[k][0]);
        if ( id1 < 0 ) {
          printf("Particle %d not in geant\n",mode[k][0]);
          return false;
        }
      }  
      else id1=0;

      if(mode[k][1]!=0) {
        id2= IdFromPDG(mode[k][1]);
        if ( id2 < 0 ) {
          printf("Particle %d not in geant\n",mode[k][1]);
          return false;
        }
      }  
      else id2=0;
      
      if(mode[k][2]!=0) {
        id3= IdFromPDG(mode[k][2]);
        if ( id3 < 0 ) {
          printf("Particle %d not in geant\n",mode[k][1]);
          return false;
        }
      }  
      else id3=0;
      g3mode[k]=id1 + id2* 100+ id3 * 10000 ;
      
   }                                      
   Gsdk(IdFromPDG(pdg), bratio, g3mode);  
   return kTRUE;
}                                

//______________________________________________________________________
void TGeant3::Vname(const char *name, char *vname)
{
  //
  //  convert name to upper case. Make vname at least 4 chars
  //
  Int_t l = strlen(name);
  Int_t i;
  l = l < 4 ? l : 4;
  for (i=0;i<l;i++) vname[i] = toupper(name[i]);
  for (i=l;i<4;i++) vname[i] = ' ';
  vname[4] = 0;
}

//______________________________________________________________________
void TGeant3::Ertrgo()
{
  //
  // Perform the tracking of the track Track parameters are in VECT
  //
  ertrgo();
}

//______________________________________________________________________
void TGeant3::Ertrak(const Float_t *x1, const Float_t *p1,
			const Float_t *x2, const Float_t *p2,
			Int_t ipa,  Option_t *chopt)
{
  //************************************************************************
  //*                                                                      *
  //*          Perform the tracking of the track from point X1 to          *
  //*                    point X2                                          *
  //*          (Before calling this routine the user should also provide   *
  //*                    the input informations in /EROPTS/ and /ERTRIO/   *
  //*                    using subroutine EUFIL(L/P/V)                     *
  //*                 X1       - Starting coordinates (Cartesian)          *
  //*                 P1       - Starting 3-momentum  (Cartesian)          *
  //*                 X2       - Final coordinates    (Cartesian)          *
  //*                 P2       - Final 3-momentum     (Cartesian)          *
  //*                 IPA      - Particle code (a la GEANT) of the track   *
  //*                                                                      *
  //*                 CHOPT                                                *
  //*                     'B'   'Backward tracking' - i.e. energy loss     *
  //*                                        added to the current energy   *
  //*                     'E'   'Exact' calculation of errors assuming     *
  //*                                        helix (i.e. path-length not   *
  //*                                        assumed as infinitesimal)     *
  //*                     'L'   Tracking up to prescribed Lengths reached  *
  //*                     'M'   'Mixed' prediction (not yet coded)         *
  //*                     'O'   Tracking 'Only' without calculating errors *
  //*                     'P'   Tracking up to prescribed Planes reached   *
  //*                     'V'   Tracking up to prescribed Volumes reached  *
  //*                     'X'   Tracking up to prescribed Point approached *
  //*                                                                      *
  //*                Interface with GEANT :                                *
  //*             Track parameters are in /CGKINE/ and /GCTRAK/            *
  //*                                                                      *
  //*          ==>Called by : USER                                         *
  //*             Authors   M.Maire, E.Nagy  ********//*                   *
  //*                                                                      *
  //************************************************************************
  ertrak(x1,p1,x2,p2,ipa,PASSCHARD(chopt) PASSCHARL(chopt));
}

void TGeant3::Erxyzc(){
//
//    ******************************************************************
//    *                                                                *
//    *        Print track and volume parameters at current point      *
//    *                                                                *
//    *    ==>Called by : <USER,EUSTEP>                                *
//    *       Author    R.Brun  *********                              *
//    *                                                                *
//    ******************************************************************
//


  erxyzc();
}



void TGeant3::Eufill(Int_t n,Float_t *ein,Float_t *xlf){

// C.    ******************************************************************
// C.    *                                                                *
// C.    *    User routine to fill the input values of the commons :      *
// C.    *               /EROPTS/, /EROPTC/ and /ERTRIO/ for CHOPT = 'L'  *
// C.    *         N     Number of predictions where to store results     *
// C.    *         EIN   Input error matrix                               *
// C.    *         XLF   Defines the tracklengths which if passed the     *
// C.    *                      result should be stored                   *
// C.    *                                                                *
// C.    *                                                                *
// C.    *    ==>Called by : USER (before calling ERTRAK)                 *
// C.    *       Author    M.Maire, E.Nagy  *********                     *
// C.    *                                                                *
// C.    ******************************************************************
   for(Int_t i=0;i<15;i++) fErtrio->errin[i]=ein[i]; 
   const Int_t mxpred=10;
   if (n<mxpred) {
      fErtrio->nepred=n;
    } else {
     fErtrio->nepred=mxpred;
   } 
   for(Int_t i=0;i<15;i++) fErtrio->errin[i]=ein[i]; 
   for(Int_t i=0;i<fErtrio->nepred;i++) fEropts->erleng[i]=xlf[i]; 
//  eufill(n,ein,xlf);
}

void TGeant3::Eufilp(const Int_t n, Float_t *ein,
			Float_t *pli, Float_t *plf)
{
  //    ******************************************************************
  //    *                                                                *
  //    *    User routine to fill the input values of the commons :      *
  //    *               /EROPTS/, /EROPTC/ and /ERTRIO/ for CHOPT = 'P'  *
  //    *         N     Number of predictions where to store results     *
  //    *         EIN   Input error matrix (in the 'Plane' system )      *
  //    *         PLI   Defines the start plane                          *
  //    *                      PLI(3,1) - and                            *
  //    *                      PLI(3,2) - 2 unit vectors in the plane    *
  //    *         PLF   Defines the end plane                            *
  //    *                      PLF(3,1,I) - and                          *
  //    *                      PLF(3,2,I) - 2 unit vectors in the plane  *
  //    *                      PLF(3,3,I) - point on the plane           *
  //    *                                   at intermediate point I      *
  //    *                                                                *
  //    *    ==>Called by : USER (before calling ERTRAK)                 *
  //    *       Author    M.Maire, E.Nagy  *********                     *
  //    *                                                                *
  //    ******************************************************************
   for(Int_t i=0;i<15;i++) fErtrio->errin[i]=ein[i]; 
   const Int_t mxpred=10;
   if (n<mxpred) {
      fErtrio->nepred=n;
    } else {
     fErtrio->nepred=mxpred;
   } 
   for(Int_t i=0;i<6;i++) fEropts->erpli[i]=pli[i]; 

   for (Int_t j=0;j<n;j++) {
     for(Int_t i=0;i<9;i++) {
       fEropts->erplo[i+12*j]=plf[i+12*j]; 
     }
     TVector3 v1(fEropts->erplo[0+12*j],fEropts->erplo[1+12*j],fEropts->erplo[2+12*j]);
     TVector3 v2(fEropts->erplo[3+12*j],fEropts->erplo[4+12*j],fEropts->erplo[5+12*j]);
     TVector3 v3=v1.Cross(v2);
     fEropts->erplo[9]=v3(0);
     fEropts->erplo[10]=v3(1);
     fEropts->erplo[11]=v3(2);
   }


}
void TGeant3::Eufilv(Int_t n, Float_t *ein,
			Char_t *namv, Int_t *numv,Int_t *iovl)
{

  //    ******************************************************************
  //    *                                                                *
  //    *    User routine to fill the input values of the commons :      *
  //    *               /EROPTS/, /EROPTC/ and /ERTRIO/ for CHOPT = 'V'  *
  //    *         N     Number of predictions where to store results     *
  //    *         EIN   Input error matrix                               *
  //    *        CNAMV  Volume name of the prediction                    *
  //    *        NUMV   Volume number (if 0 = all volumes)               *
  //    *        IOVL   = 1  prediction when entering in the volume      *
  //    *               = 2  prediction when leaving the volume          *
  //    *                                                                *
  //    *    ==>Called by : USER (before calling ERTRAK)                 *
  //    *       Author    M.Maire, E.Nagy  *********                     *
  //    *                                                                *
  //    ******************************************************************

  for(Int_t i=0;i<15;i++) fErtrio->errin[i]=ein[i]; 
   const Int_t mxpred=15;
   if (n<mxpred) {
      fErtrio->nepred=n;
    } else {
     fErtrio->nepred=mxpred;
   } 
   
   for(Int_t i=0;i<fErtrio->nepred;i++) {
     fEropts->nameer[i]=*((int*)namv);
     fEropts->iovler[i]=iovl[i];
     fEropts->numver[i]=numv[i];
   }
}
//______________________________________________________________________
void TGeant3::Trscsd(Float_t *pc,Float_t *rc,Float_t *pd,Float_t *rd,Float_t *h,Float_t ch,Int_t ierr,Float_t spu,Float_t *dj,Float_t *dk){

//       SUBROUTINE TRSCSD(PC,RC,PD,RD,H,CH,IERR,SPU,DJ,DK)
// ******************************************************************
//  *** TRANSFORMS ERROR MATRIX
//    FROM   SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
//       TO         VARIABLES (1/P,V',W',V,W)
// 
//      Authors: A. Haas and W. Wittek
//  *** PC(3)     1/P,LAMBDA,PHI                          INPUT
//      PD(3)     1/P,V',W'                              OUTPUT
//      H(3)      MAGNETIC FIELD                          INPUT
//      RC(15)    ERROR MATRIX IN   SC   VARIABLES        INPUT     (TRIANGLE)
//      RD(15)    ERROR MATRIX IN 1/P,V',W',V,W          OUTPUT     (TRIANGLE)
//      CH        CHARGE OF PARTICLE                      INPUT
//                CHARGE AND MAGNETIC FIELD ARE NEEDED
//                FOR CORRELATION TERMS (V',YT),(V',ZT),(W',YT),(W',ZT)
//                THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
//                TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
//                AND RD FOR FIXED U
//      DJ(3)     UNIT VECTOR IN V-DIRECTION
//      DK(3)     UNIT VECTOR IN W-DIRECTION    OF DETECTOR SYSTEM
// 
//      IERR  =   1       PARTICLE MOVES PERPENDICULAR TO U-AXIS
//                       ( V',W' ARE NOT DEFINED )
//      SPU       SIGN OF U-COMPONENT OF PARTICLE MOMENTUM   OUTPUT
// ******************************************************************
  printf("%d\n",ierr);
  trscsd(pc,rc,pd,rd,h,ch,ierr,spu,dj,dk);
}
//______________________________________________________________________
void TGeant3::Trsdsc(Float_t *pd,Float_t *rd,Float_t *pc,Float_t *rc,Float_t *h,Float_t *ch,Int_t *ierr,Float_t *spu,Float_t *dj,Float_t *dk) {
// ******************************************************************
//       SUBROUTINE TRSDSC(PD,RD,PC,RC,H,CH,IERR,SPU,DJ,DK)
// 
//  *** TRANSFORMS ERROR MATRIX
//      FROM        VARIABLES (1/P,V',W',V,W)
//       TO    SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
//      Authors: A. Haas and W. Wittek
//  *** PD(3)     1/P,V',W'                               INPUT
//      PC(3)     1/P,LAMBDA,PHI                         OUTPUT
//      H(3)      MAGNETIC FIELD                          INPUT
//      RD(15)    ERROR MATRIX IN 1/P,V',W',V,W           INPUT      (TRIANGLE)
//      RC(15)    ERROR MATRIX IN   SC   VARIABLES       OUTPUT      (TRIANGLE)
//      CH        CHARGE OF PARTICLE                      INPUT
//                CHARGE AND MAGNETIC FIELD ARE NEEDED
//                FOR CORRELATION TERMS (LAMBDA,V),(LAMBDA,W),(PHI,V),(PHI,W)
//                THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
//                TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
//                AND RD FOR FIXED U
//      DJ(3)     UNIT VECTOR IN V-DIRECTION
//      DK(3)     UNIT VECTOR IN W-DIRECTION    OF DETECTOR SYSTEM
// 
//      IERR              NOT USED
//      SPU       SIGN OF U-COMPONENT OF PARTICLE MOMENTUM    INPUT
// ******************************************************************
 trsdsc(pd,rd,pc,rc,h,ch,ierr,spu,dj,dk);
}
//______________________________________________________________________
void TGeant3::Trscsp(Float_t *pc,Float_t *rc,Float_t *ps,Float_t *rs,Float_t *h,Float_t *ch,Int_t *ierr, Float_t *spx){
// ******************************************************************
//       SUBROUTINE TRSCSP(PC,RC,PS,RS,H,CH,IERR,SPX)
// 
//  *** TRANSFORMS ERROR MATRIX
//      FROM   SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
//       TO  SPLINE VARIABLES (1/P,Y',Z',Y,Z)
// 
//      Authors: A. Haas and W. Wittek
// 
// 
//  *** PC(3)     1/P,LAMBDA,PHI                          INPUT
//      PS(3)     1/P,Y',Z'                              OUTPUT
//      H(3)      MAGNETIC FIELD                          INPUT
//      RC(15)    ERROR MATRIX IN   SC   VARIABLES        INPUT     (TRIANGLE)
//      RS(15)    ERROR MATRIX IN SPLINE VARIABLES       OUTPUT     (TRIANGLE)
//      CH        CHARGE OF PARTICLE                      INPUT
//                CHARGE AND MAGNETIC FIELD ARE NEEDED
//                FOR CORRELATION TERMS (Y',YT),(Y',ZT),(Z',YT),(Z',ZT)
//                THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
//                TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
//                AND RS FOR FIXED X
// 
//      IERR  =   1       PARTICLE MOVES PERPENDICULAR TO X-AXIS
//                       ( Y',Z' ARE NOT DEFINED )
//      SPX       SIGN OF X-COMPONENT OF PARTICLE MOMENTUM   OUTPUT
// ******************************************************************
  trscsp(pc,rc,ps,rs,h,ch,ierr,spx);
}
//______________________________________________________________________
void TGeant3::Trspsc(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,Float_t *ch,Int_t *ierr,Float_t *spx) {

//     ******************************************************************
//       SUBROUTINE TRSPSC(PS,RS,PC,RC,H,CH,IERR,SPX)
// 
//  *** TRANSFORMS ERROR MATRIX
//      FROM SPLINE VARIABLES (1/P,Y',Z',Y,Z)
//       TO    SC   VARIABLES (1/P,LAMBDA,PHI,YT,ZT)
// 
//      Authors: A. Haas and W. Wittek
// 
// 
//  *** PS(3)     1/P,Y',Z'                               INPUT
//      PC(3)     1/P,LAMBDA,PHI                         OUTPUT
//      H(3)      MAGNETIC FIELD                          INPUT
//      RS(15)    ERROR MATRIX IN SPLINE VARIABLES        INPUT      (TRIANGLE)
//      RC(15)    ERROR MATRIX IN   SC   VARIABLES       OUTPUT      (TRIANGLE)
//      CH        CHARGE OF PARTICLE                      INPUT
//                CHARGE AND MAGNETIC FIELD ARE NEEDED
//                FOR CORRELATION TERMS (LAMBDA,Y),(LAMBDA,Z),(PHI,Y),(PHI,Z)
//                THESE CORRELATION TERMS APPEAR BECAUSE RC IS ASSUMED
//                TO BE THE ERROR MATRIX FOR FIXED S (PATH LENGTH)
//                AND RS FOR FIXED X
// 
//      IERR              NOT USED
//      SPX       SIGN OF X-COMPONENT OF PARTICLE MOMENTUM    INPUT
// 
//     ******************************************************************

 trspsc(ps,rs,pc,rc,h,ch,ierr,spx);

}


//______________________________________________________________________
void TGeant3::WriteEuclid(const char* filnam, const char* topvol,
			  Int_t number, Int_t nlevel)
{
  //
  //
  //     ******************************************************************
  //     *                                                                *
  //     *  Write out the geometry of the detector in EUCLID file format  *
  //     *                                                                *
  //     *       filnam : will be with the extension .euc                 *
  //     *       topvol : volume name of the starting node                *
  //     *       number : copy number of topvol (relevant for gsposp)     *
  //     *       nlevel : number of  levels in the tree structure         *
  //     *                to be written out, starting from topvol         *
  //     *                                                                *
  //     *       Author : M. Maire                                        *
  //     *                                                                *
  //     ******************************************************************
  //
  //     File filnam.tme is written out with the definitions of tracking
  //     medias and materials.
  //     As to restore original numbers for materials and medias, program
  //     searches in the file euc_medi.dat and comparing main parameters of
  //     the mat. defined inside geant and the one in file recognizes them
  //     and is able to take number from file. If for any material or medium,
  //     this procedure fails, ordering starts from 1.
  //     Arrays IOTMED and IOMATE are used for this procedure
  //
  const char kShape[][5]={"BOX ","TRD1","TRD2","TRAP","TUBE","TUBS","CONE",
			 "CONS","SPHE","PARA","PGON","PCON","ELTU","HYPE",
			 "GTRA","CTUB"};
  Int_t i, end, itm, irm, jrm, k, nmed;
  Int_t imxtmed=0;
  Int_t imxmate=0;
  FILE *lun;
  char *filext, *filetme;
  char natmed[21], namate[21];
  char natmedc[21], namatec[21];
  char key[5], name[5], mother[5], konly[5];
  char card[133];
  Int_t iadvol, iadtmd, iadrot, nwtot, iret;
  Int_t mlevel, numbr, natt, numed, nin, ndata;
  Int_t iname, ivo, ish, jvo, nvstak, ivstak;
  Int_t jdiv, ivin, in, jin, jvin, irot;
  Int_t jtm, imat, jma, flag=0, imatc;
  Float_t az, dens, radl, absl, a, step, x, y, z;
  Int_t npar, ndvmx, left;
  Float_t zc, densc, radlc, abslc, c0, tmaxfd;
  Int_t nparc, numb;
  Int_t iomate[100], iotmed[100];
  Float_t par[100], att[20], ubuf[50];
  Float_t *qws;
  Int_t   *iws;
  Int_t level, ndiv, iaxe;
  Int_t itmedc, nmatc, isvolc, ifieldc, nwbufc, isvol, nmat, ifield, nwbuf;
  Float_t fieldmc, tmaxfdc, stemaxc, deemaxc, epsilc, stminc, fieldm;
  Float_t tmaxf, stemax, deemax, epsil, stmin;
  const char *k10000="!\n%s\n!\n";
  //Open the input file
  end=strlen(filnam);
  for(i=0;i<end;i++) if(filnam[i]=='.') {
    end=i;
    break;
  }
  filext=new char[end+5];
  filetme=new char[end+5];
  strncpy(filext,filnam,end);
  strncpy(filetme,filnam,end);
  //
  // *** The output filnam name will be with extension '.euc'
  strcpy(&filext[end],".euc");
  strcpy(&filetme[end],".tme");
  lun=fopen(filext,"w");
  //
  // *** Initialization of the working space
  iadvol=fGcnum->nvolum;
  iadtmd=iadvol+fGcnum->nvolum;
  iadrot=iadtmd+fGcnum->ntmed;
  if(fGclink->jrotm) {
    fGcnum->nrotm=fZiq[fGclink->jrotm-2];
  } else {
    fGcnum->nrotm=0;
  }
  nwtot=iadrot+fGcnum->nrotm;
  qws = new float[nwtot+1];
  for (i=0;i<nwtot+1;i++) qws[i]=0;
  iws = (Int_t*) qws;
  mlevel=nlevel;
  if(nlevel==0) mlevel=20;
  //
  // *** find the top volume and put it in the stack
  numbr = number>0 ? number : 1;
  Gfpara(topvol,numbr,1,npar,natt,par,att);
  if(npar <= 0) {
    printf(" *** GWEUCL *** top volume : %s number : %3d can not be "
           "a valid root\n", topvol, numbr);
    return;
  }
  //
  // ***  authorized shape ?
  strncpy((char *)&iname, topvol, 4);
  ivo=0;
  for(i=1; i<=fGcnum->nvolum; i++) if(fZiq[fGclink->jvolum+i]==iname) {
    ivo=i;
    break;
  }
  jvo = fZlq[fGclink->jvolum-ivo];
  ish = Int_t (fZq[jvo+2]);
  if(ish > 12) {
    printf(" *** GWEUCL *** top volume : %s number : %3d can not be "
           "a valid root\n",topvol, numbr);
  }
  //
  level = 1;
  nvstak = 1;
  iws[nvstak]     = ivo;
  iws[iadvol+ivo] = level;
  ivstak = 0;
  //
  //*** flag all volumes and fill the stack
  //
 L10:
  //
  //    pick the next volume in stack
  ivstak += 1;
  ivo   = TMath::Abs(iws[ivstak]);
  jvo   = fZlq[fGclink->jvolum - ivo];
  //
  //     flag the tracking medium
  numed =  Int_t (fZq[jvo + 4]);
  iws[iadtmd + numed] = 1;
  //
  //    get the daughters ...
  level = iws[iadvol+ivo];
  if (level < mlevel) {
    level +=  1;
    nin = Int_t (fZq[jvo + 3]);
    //
    //       from division ...
    if (nin < 0) {
      jdiv = fZlq[jvo  - 1];
      ivin =  Int_t (fZq[jdiv + 2]);
      nvstak += 1;
      iws[nvstak]      = -ivin;
      iws[iadvol+ivin] =  level;
      //
      //       from position ...
    } else if (nin > 0) {
      for(in=1; in<=nin; in++) {
	jin  = fZlq[jvo - in];
	ivin =  Int_t (fZq[jin + 2 ]);
	jvin = fZlq[fGclink->jvolum - ivin];
	ish  =  Int_t (fZq[jvin + 2]);
	//              authorized shape ?
	if (ish <= 12) {
	  //                 not yet flagged ?
	  if (iws[iadvol+ivin]==0) {
	    nvstak += 1;
	    iws[nvstak]      = ivin;
	    iws[iadvol+ivin] = level;
	  }
	  //                 flag the rotation matrix
	  irot = Int_t ( fZq[jin + 4 ]);
	  if (irot > 0) iws[iadrot+irot] = 1;
	}
      }
    }
  }
  //
  //     next volume in stack ?
  if (ivstak < nvstak) goto L10;
  //
  // *** restore original material and media numbers
  // file euc_medi.dat is needed to compare materials and medias
  //
  FILE* luncor=fopen("euc_medi.dat","r");
  //
  if(luncor) {
    for(itm=1; itm<=fGcnum->ntmed; itm++) {
      if (iws[iadtmd+itm] > 0) {
	jtm = fZlq[fGclink->jtmed-itm];
	strncpy(natmed,(char *)&fZiq[jtm+1],20);
	imat =  Int_t (fZq[jtm+6]);
	jma  = fZlq[fGclink->jmate-imat];
	if (jma <= 0) {
	  printf(" *** GWEUCL *** material not defined for tracking medium "
              "%5i %s\n",itm,natmed);
	  flag=1;
	} else {
	  strncpy(namate,(char *)&fZiq[jma+1],20);
	}
	//*
	//** find the material original number
	rewind(luncor);
      L23:
	iret=fscanf(luncor,"%4s,%130s",key,card);
	if(iret<=0) goto L26;
	flag=0;
	if(!strcmp(key,"MATE")) {
	  sscanf(card,"%d %s %f %f %f %f %f %d",&imatc,namatec,&az,&zc,
              &densc,&radlc,&abslc,&nparc);
	  Gfmate(imat,namate,a,z,dens,radl,absl,par,npar);
	  if(!strcmp(namatec,namate)) {
	    if(az==a && zc==z && densc==dens && radlc==radl
	       && abslc==absl && nparc==nparc) {
	      iomate[imat]=imatc;
	      flag=1;
	      printf("*** GWEUCL *** material : %3d '%s' restored as %3d\n",
                  imat,namate,imatc);
	    } else {
	      printf("*** GWEUCL *** different definitions for material: %s\n",
                  namate);
	    }
	  }
	}
	if(strcmp(key,"END") && !flag) goto L23;
	if (!flag) {
	  printf("*** GWEUCL *** cannot restore original number for "
              "material: %s\n",namate);
	}
	//*
	//*
	//***  restore original tracking medium number
	rewind(luncor);
      L24:
	iret=fscanf(luncor,"%4s,%130s",key,card);
	if(iret<=0) goto L26;
	flag=0;
	if (!strcmp(key,"TMED")) {
	  sscanf(card,"%d %s %d %d %d %f %f %f %f %f %f %d\n",
		 &itmedc,natmedc,&nmatc,&isvolc,&ifieldc,&fieldmc,
		 &tmaxfdc,&stemaxc,&deemaxc,&epsilc,&stminc,&nwbufc);
	  Gftmed(itm,natmed,nmat,isvol,ifield,fieldm,tmaxf,stemax,deemax,
			epsil,stmin,ubuf,&nwbuf);
	  if(!strcmp(natmedc,natmed)) {
	    if (iomate[nmat]==nmatc && nwbuf==nwbufc) {
	      iotmed[itm]=itmedc;
	      flag=1;
	      printf("*** GWEUCL *** medium   : %3d '%20s' restored as %3d\n",
		     itm,natmed,itmedc);
	    } else {
	      printf("*** GWEUCL *** different definitions for tracking "
                  "medium: %s\n",natmed);
	    }
	  }
	}
	if(strcmp(key,"END") && !flag) goto L24;
	if(!flag) {
	  printf("cannot restore original number for medium : %s\n",natmed);
	  goto L27;
	}
      }
    }
    goto L29;
    //*
  }
 L26:   printf("*** GWEUCL *** cannot read the data file\n");
 L27:   flag=2;
 L29:   if(luncor) fclose (luncor);
  //
  //
  // *** write down the tracking medium definition
  //
  strcpy(card,"!       Tracking medium");
  fprintf(lun,k10000,card);
  //
  for(itm=1;itm<=fGcnum->ntmed;itm++) {
    if (iws[iadtmd+itm]>0) {
      jtm  = fZlq[fGclink->jtmed-itm];
      strncpy(natmed,(char *)&fZiq[jtm+1],20);
      natmed[20]='\0';
      imat =  Int_t (fZq[jtm+6]);
      jma  = fZlq[fGclink->jmate-imat];
      //*  order media from one, if comparing with database failed
      if (flag==2) {
	iotmed[itm]=++imxtmed;
	iomate[imat]=++imxmate;
      }
      //*
      if(jma<=0) {
	strcpy(namate,"                  ");
	printf(" *** GWEUCL *** material not defined for tracking "
            "medium %5d %s\n", itm,natmed);
      } else {
	strncpy(namate,(char *)&fZiq[jma+1],20);
	namate[20]='\0';
      }
      fprintf(lun,"TMED %3d '%20s' %3d '%20s'\n",iotmed[itm],natmed,
              iomate[imat],namate);
    }
  }
  //*
      //* *** write down the rotation matrix
  //*
  strcpy(card,"!       Reperes");
  fprintf(lun,k10000,card);
  //
  for(irm=1;irm<=fGcnum->nrotm;irm++) {
    if (iws[iadrot+irm]>0) {
      jrm  = fZlq[fGclink->jrotm-irm];
      fprintf(lun,"ROTM %3d",irm);
      for(k=11;k<=16;k++) fprintf(lun," %8.3f",fZq[jrm+k]);
      fprintf(lun,"\n");
    }
  }
  //*
  //* *** write down the volume definition
  //*
  strcpy(card,"!       Volumes");
  fprintf(lun,k10000,card);
  //*
  for(ivstak=1;ivstak<=nvstak;ivstak++) {
    ivo = iws[ivstak];
    if (ivo>0) {
      strncpy(name,(char *)&fZiq[fGclink->jvolum+ivo],4);
      name[4]='\0';
      jvo  = fZlq[fGclink->jvolum-ivo];
      ish   = Int_t (fZq[jvo+2]);
      nmed  = Int_t (fZq[jvo+4]);
      npar  = Int_t (fZq[jvo+5]);
      if (npar>0) {
	if (ivstak>1) for(i=0;i<npar;i++) par[i]=fZq[jvo+7+i];
	Gckpar (ish,npar,par);
	fprintf(lun,"VOLU '%4s' '%4s' %3d %3d\n",name,kShape[ish-1],
             iotmed[nmed],npar);
	for(i=0;i<(npar-1)/6+1;i++) {
	  fprintf(lun,"     ");
	  left=npar-i*6;
	  for(k=0;k<(left<6?left:6);k++) fprintf(lun," %11.5f",par[i*6+k]);
	  fprintf(lun,"\n");
	}
      } else {
	fprintf(lun,"VOLU '%4s' '%4s' %3d %3d\n",name,kShape[ish-1],
             iotmed[nmed],npar);
      }
    }
  }
  //*
  //* *** write down the division of volumes
  //*
  fprintf(lun,k10000,"!       Divisions");
  for(ivstak=1;ivstak<=nvstak;ivstak++) {
    ivo = TMath::Abs(iws[ivstak]);
    jvo  = fZlq[fGclink->jvolum-ivo];
    ish  = Int_t (fZq[jvo+2]);
    nin  = Int_t (fZq[jvo+3]);
    //*        this volume is divided ...
    if (nin<0) {
      jdiv = fZlq[jvo-1];
      iaxe = Int_t ( fZq[jdiv+1]);
      ivin = Int_t ( fZq[jdiv+2]);
      ndiv = Int_t ( fZq[jdiv+3]);
      c0   =  fZq[jdiv+4];
      step =  fZq[jdiv+5];
      jvin = fZlq[fGclink->jvolum-ivin];
      nmed = Int_t ( fZq[jvin+4]);
      strncpy(mother,(char *)&fZiq[fGclink->jvolum+ivo ],4);
      mother[4]='\0';
      strncpy(name,(char *)&fZiq[fGclink->jvolum+ivin],4);
      name[4]='\0';
      if ((step<=0.)||(ish>=11)) {
	//*              volume with negative parameter or gsposp or pgon ...
	fprintf(lun,"DIVN '%4s' '%4s' %3d %3d\n",name,mother,ndiv,iaxe);
      } else if ((ndiv<=0)||(ish==10)) {
	//*              volume with negative parameter or gsposp or para ...
	ndvmx = TMath::Abs(ndiv);
	fprintf(lun,"DIVT '%4s' '%4s' %11.5f %3d %3d %3d\n",
		name,mother,step,iaxe,iotmed[nmed],ndvmx);
      } else {
	//*              normal volume : all kind of division are equivalent
	fprintf(lun,"DVT2 '%4s' '%4s' %11.5f %3d %11.5f %3d %3d\n",
		name,mother,step,iaxe,c0,iotmed[nmed],ndiv);
      }
    }
  }
  //*
  //* *** write down the the positionnement of volumes
  //*
  fprintf(lun,k10000,"!       Positionnements\n");
  //
  for(ivstak = 1;ivstak<=nvstak;ivstak++) {
    ivo = TMath::Abs(iws[ivstak]);
    strncpy(mother,(char*)&fZiq[fGclink->jvolum+ivo ],4);
    mother[4]='\0';
    jvo  = fZlq[fGclink->jvolum-ivo];
    nin  = Int_t( fZq[jvo+3]);
    //*        this volume has daughters ...
    if (nin>0) {
      for (in=1;in<=nin;in++) {
	jin  = fZlq[jvo-in];
	ivin =  Int_t (fZq[jin +2]);
	numb =  Int_t (fZq[jin +3]);
	irot =  Int_t (fZq[jin +4]);
	x    =  fZq[jin +5];
	y    =  fZq[jin +6];
	z    =  fZq[jin +7];
	strcpy(konly,"ONLY");
	if (fZq[jin+8]!=1.) strcpy(konly,"MANY");
	strncpy(name,(char*)&fZiq[fGclink->jvolum+ivin],4);
	name[4]='\0';
	jvin = fZlq[fGclink->jvolum-ivin];
	ish  = Int_t (fZq[jvin+2]);
	//*              gspos or gsposp ?
	ndata = fZiq[jin-1];
	if (ndata==8) {
	  fprintf(lun,"POSI '%4s' %4d '%4s' %11.5f %11.5f %11.5f %3d '%4s'\n",
		  name,numb,mother,x,y,z,irot,konly);
	} else {
	  npar =  Int_t (fZq[jin+9]);
	  for(i=0;i<npar;i++) par[i]=fZq[jin+10+i];
	  Gckpar (ish,npar,par);
	  fprintf(lun,"POSP '%4s' %4d '%4s' %11.5f %11.5f %11.5f %3d '%4s' %3d\n",
		  name,numb,mother,x,y,z,irot,konly,npar);
	  fprintf(lun,"     ");
	  for(i=0;i<npar;i++) fprintf(lun," %11.5f",par[i]);
	  fprintf(lun,"\n");
	}
      }
    }
  }
  //*
  fprintf(lun,"END\n");
  fclose(lun);
  //*
  //****** write down the materials and medias *****
  //*
  lun=fopen(filetme,"w");
  //*
  for(itm=1;itm<=fGcnum->ntmed;itm++) {
    if (iws[iadtmd+itm]>0) {
      jtm  = fZlq[fGclink->jtmed-itm];
      strncpy(natmed,(char*)&fZiq[jtm+1],4);
      imat =  Int_t (fZq[jtm+6]);
      jma  =  Int_t (fZlq[fGclink->jmate-imat]);
      //*  material
      Gfmate (imat,namate,a,z,dens,radl,absl,par,npar);
      fprintf(lun,"MATE %4d '%20s'%11.5E %11.5E %11.5E %11.5E %11.5E %3d\n",
	     iomate[imat],namate,a,z,dens,radl,absl,npar);
      //*
      if (npar>0) {
	  fprintf(lun,"     ");
	  for(i=0;i<npar;i++) fprintf(lun," %11.5f",par[i]);
	  fprintf(lun,"\n");
      }
      //*  medium
      Gftmed(itm,natmed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,deemax,
             epsil,stmin,par,&npar);
      fprintf(lun,"TMED %4d '%20s' %3d %1d %3d %11.5f %11.5f %11.5f "
              "%11.5f %11.5f %11.5f %3d\n",
              iotmed[itm],natmed,iomate[nmat],isvol,ifield,
              fieldm,tmaxfd,stemax,deemax,epsil,stmin,npar);
      //*
      if (npar>0) {
	  fprintf(lun,"     ");
	  for(i=0;i<npar;i++) fprintf(lun," %11.5f",par[i]);
	  fprintf(lun,"\n");
      }

    }
  }
  fprintf(lun,"END\n");
  fclose(lun);
  printf(" *** GWEUCL *** file: %s is now written out\n",filext);
  printf(" *** GWEUCL *** file: %s is now written out\n",filetme);
  // Clean up
  delete [] filext;
  delete [] filetme;
  delete [] qws;
  iws=0;
  return;
}

//______________________________________________________________________
Int_t  TGeant3::TransportMethod(TMCParticleType particleType) const
{
//
// Returns G3 transport method code for the specified MCParticleType
// ---

  switch (particleType) {
    case kPTGamma:    return 1;
    case kPTElectron: return 2;
    case kPTNeutron:  return 3;
    case kPTHadron:   return 4;
    case kPTMuon:     return 5;
    case kPTGeantino: return 6;
    case kPTOpticalPhoton: return 7;
    case kPTIon:      return 8;
    default:          return -1;
  }
}

//______________________________________________________________________
TMCParticleType  TGeant3::ParticleType(Int_t itrtyp) const
{
//
// Returns MCParticleType for the specified G3 transport method code
// ---

  switch (itrtyp) {
    case 1:  return  kPTGamma;
    case 2:  return  kPTElectron;
    case 3:  return  kPTNeutron;
    case 4:  return  kPTHadron;
    case 5:  return  kPTMuon;
    case 6:  return  kPTGeantino;
    case 7:  return  kPTOpticalPhoton;
    case 8:  return  kPTIon;
    default: return  kPTUndefined;
  }
}

//______________________________________________________________________
TString  TGeant3::ParticleClass(TMCParticleType particleType) const
{
//
// Returns particle class name (used in TDatabasePDG) for
// the specified MCParticleType
// ---

  // CHECK
  switch (particleType) {
    case kPTGamma:    return TString("Photon");
    case kPTElectron: return TString("Lepton");
    case kPTNeutron:  return TString("Hadron");
    case kPTHadron:   return TString("Hadron");
    case kPTMuon:     return TString("Lepton");
    case kPTGeantino: return TString("Special");
    case kPTIon:      return TString("Ion");
    case kPTOpticalPhoton: return TString("Photon");
    default:          return TString("Unknown");
  }
}

//______________________________________________________________________
void TGeant3::FinishGeometry()
{
  //
  // Finalize geometry construction
  //

  //Close the geometry structure
  if (gDebug > 0) printf("FinishGeometry, calling ggclos\n");
  Ggclos();


  //  gROOT->GetListOfBrowsables()->Add(gGeoManager);
  if (gDebug > 0) printf("FinishGeometry, calling SetColors\n");

  //Create the color table
  SetColors();
  if (gDebug > 0) printf("FinishGeometry, returning\n");
}

//______________________________________________________________________
void TGeant3::Init()
{
    //
    //=================Create Materials and geometry
    //

    //  Some default settings, if not changed by user
    if (!TestBit(kTRIG)) SetTRIG(1);       // Number of events to be processed
    if (!TestBit(kSWIT)) SetSWIT(4, 10);   //
    if (!TestBit(kDEBU)) SetDEBU(0, 0, 1); //
    if (!TestBit(kAUTO)) SetAUTO(1);       // Select automatic STMIN etc... 
                                           // calc. (AUTO 1) or manual (AUTO 0)
    if (!TestBit(kABAN)) SetABAN(0);       // Restore 3.16 behaviour for 
                                           // abandoned tracks
    if (!TestBit(kOPTI)) SetOPTI(2);       // Select optimisation level for 
                                           // GEANT geometry searches (0,1,2)
    if (!TestBit(kERAN)) SetERAN(5.e-7);   //

    DefineParticles();
    fApplication->AddParticles();
    fApplication->AddIons();
    fApplication->ConstructGeometry();
    FinishGeometry();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,01,1)
    fApplication->ConstructOpGeometry();
#endif
    fApplication->InitGeometry();
}

//____________________________________________________________________________
Bool_t TGeant3::ProcessRun(Int_t nevent)
{
  //
  // Process the run and return true if run has finished successfully,
  // return false in other cases (run aborted by user)

  Int_t todo = TMath::Abs(nevent);
  for (Int_t i=0; i<todo; i++) {
     // Process one run (one run = one event)
     fGcflag->idevt  = i;
     fGcflag->ievent = i+1;
     if (fStopRun) break;
     fApplication->BeginEvent();
     if (fStopRun) break;
     ProcessEvent();
     if (fStopRun) break;
     fApplication->FinishEvent();
     if (fStopRun) break;
  }

  if (fStopRun) printf(" **** Run stopped ***\n");

  Bool_t returnValue = !fStopRun;
  fStopRun = kFALSE;
#ifdef STATISTICS
  printf("count_gmedia= %8d\n",count_gmedia);
  printf("count_gtmedi= %8d\n",count_gtmedi);
  printf("count_ginvol= %8d\n",count_ginvol);
  printf("count_gtnext= %8d\n",count_gtnext);
  stattree->AutoSave();
  statfile->Close();
  printf("Statistics tree saved.\n");
#endif
  return returnValue;
}

//______________________________________________________________________
void TGeant3::ProcessEvent()
{
  //
  // Process one event
  //
  Gtrigi();
  Gtrigc();
  Gtrig();
}

//______________________________________________________________________
void TGeant3::SetColors()
{
  //
  // Set the colors for all the volumes
  // this is done sequentially for all volumes
  // based on the number of their medium
  //

  Int_t kv, icol;
  Int_t jvolum=fGclink->jvolum;
  //Int_t jtmed=fGclink->jtmed;
  //Int_t jmate=fGclink->jmate;
  Int_t nvolum=fGcnum->nvolum;
  char name[5];
  //
  //    Now for all the volumes
  for(kv=1;kv<=nvolum;kv++) {
    //     Get the tracking medium
    Int_t itm=Int_t (fZq[fZlq[jvolum-kv]+4]);
    //     Get the material
    //Int_t ima=Int_t (fZq[fZlq[jtmed-itm]+6]);
    //     Get z
    //Float_t z=fZq[fZlq[jmate-ima]+7];
    //     Find color number
    //icol = Int_t(z)%6+2;
    //icol = 17+Int_t(z*150./92.);
    //icol = kv%6+2;
    icol = itm%6+2;
    strncpy(name,(char*)&fZiq[jvolum+kv],4);
    name[4]='\0';
    Gsatt(name,"COLO",icol);
  }
}

//______________________________________________________________________
void TGeant3::SetTrack(Int_t done, Int_t parent, Int_t pdg, Float_t *pmom,
		        Float_t *vpos, Float_t *polar, Float_t tof,
		        TMCProcess mech, Int_t &ntr, Float_t weight, Int_t is)
{
  //
  // Load a track on the stack
  //
  // done     0 if the track has to be transported
  //          1 if not
  // parent   identifier of the parent track. -1 for a primary
  // pdg    particle code
  // pmom     momentum GeV/c
  // vpos     position
  // polar    polarization
  // tof      time of flight in seconds
  // mecha    production mechanism
  // ntr      on output the number of the track stored
  //

  //  const Float_t tlife=0;

  //
  // Here we get the static mass
  // For MC is ok, but a more sophisticated method could be necessary
  // if the calculated mass is required
  // also, this method is potentially dangerous if the mass
  // used in the MC is not the same of the PDG database
  //
  Float_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  Float_t e=TMath::Sqrt(mass*mass+pmom[0]*pmom[0]+
			pmom[1]*pmom[1]+pmom[2]*pmom[2]);

//    printf("Loading  mass %f ene %f No %d ip %d parent %d done %d "
//           "pos %f %f %f mom %f %f %f kS %d m \n",
//	        mass,e,fNtrack,pdg,parent,done,vpos[0],vpos[1],vpos[2],
//           pmom[0],pmom[1],pmom[2],kS);


  GetStack()->PushTrack(done, parent, pdg, pmom[0], pmom[1], pmom[2], e,
                       vpos[0],vpos[1],vpos[2],tof,polar[0],polar[1],polar[2],
                       mech, ntr, weight, is);
}


//______________________________________________________________________
Float_t* TGeant3::CreateFloatArray(Float_t* array, Int_t size) const
{
// Converts Double_t* array to Float_t*,
// !! The new array has to be deleted by user.
// ---

  Float_t* floatArray;
  if (size>0) {
    floatArray = new Float_t[size];
    for (Int_t i=0; i<size; i++)
      if (array[i] >= FLT_MAX ) 
        floatArray[i] = FLT_MAX/100.;
      else	
        floatArray[i] = array[i];
  }
  else {
    //floatArray = 0;
    floatArray = new Float_t[1];
  }
  return floatArray;
}


//______________________________________________________________________
Float_t* TGeant3::CreateFloatArray(Double_t* array, Int_t size) const
{
// Converts Double_t* array to Float_t*,
// !! The new array has to be deleted by user.
// ---

  Float_t* floatArray;
  if (size>0) {
    floatArray = new Float_t[size];
    for (Int_t i=0; i<size; i++)
      if (array[i] >= FLT_MAX ) 
        floatArray[i] = FLT_MAX/100.;
      else	
        floatArray[i] = array[i];
  }
  else {
    //floatArray = 0;
    floatArray = new Float_t[1];
  }
  return floatArray;
}


//______________________________________________________________________
void TGeant3::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class TGeant3.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TVirtualMC::Streamer(R__b);
    R__b >> fNextVol;
    R__b >> fNPDGCodes;
    //R__b.ReadStaticArray(fPDGCode);
    fPDGCode.Streamer(R__b);
  } else {
    R__b.WriteVersion(TGeant3::IsA());
    TVirtualMC::Streamer(R__b);
    R__b << fNextVol;
    R__b << fNPDGCodes;
    //R__b.WriteArray(fPDGCode, fNPDGCodes);
    fPDGCode.Streamer(R__b);
  }
}

//______________________________________________________________________
//
//                 Interfaces to Fortran
//
//______________________________________________________________________


//______________________________________________________________________
extern "C" void type_of_call  rxgtrak(Int_t &mtrack,Int_t &ipart,Float_t *pmom,
                                      Float_t &e,Float_t *vpos,Float_t *polar,
                                      Float_t &tof)
{
  //
  //     Fetches next track from the ROOT stack for transport. Called by the
  //     modified version of GTREVE.
  //
  //              Track number in the ROOT stack. If MTRACK=0 no
  //      mtrack  more tracks are left in the stack to be
  //              transported.
  //      ipart   Particle code in the GEANT conventions.
  //      pmom[3] Particle momentum in GeV/c
  //      e       Particle energy in GeV
  //      vpos[3] Particle position
  //      tof     Particle time of flight in seconds
  //

  TParticle* track = gMC->GetStack()->PopNextTrack(mtrack);

  if (track) {
    // fill G3 arrays
    pmom[0] = track->Px();
    pmom[1] = track->Py();
    pmom[2] = track->Pz();
    e = track->Energy();
    vpos[0] = track->Vx();;
    vpos[1] = track->Vy();
    vpos[2] = track->Vz();
    tof = track->T();
    TVector3 pol;
    track->GetPolarisation(pol);
    polar[0] = pol.X();
    polar[1] = pol.Y();
    polar[2] = pol.Z();
    ipart = gMC->IdFromPDG(track->GetPdgCode());
  }

  mtrack++;
}


//______________________________________________________________________
extern "C" void type_of_call  rxouth ()
{
  //
  // Called by Gtreve at the end of each primary track
  //
  TVirtualMCApplication::Instance()->FinishPrimary();
}

//______________________________________________________________________
extern "C" void type_of_call  rxinh ()
{
  //
  // Called by Gtreve at the beginning of each primary track
  //
  TVirtualMCApplication::Instance()->BeginPrimary();
}

//______________________________________________________________________
void ginvol(Float_t *x, Int_t &isame)
{
   fginvol(x,isame);
}


//______________________________________________________________________
void gtmedi(Float_t *x, Int_t &numed)
{
   fgtmedi(x,numed);
#ifdef STATISTICS
   statcode = 2;
   statsame = gcchan->lsamvl;
   for (int j=0;j<6;j++) if (j <3) oldvect[j] = x[j]; else oldvect[j]=0;
   oldsafety = gctrak->safety;
   oldstep   = gctrak->step;
   sprintf(statpath,"%s",geant3->GetPath());
   statsnext=gctrak->snext;
   statsafety=gctrak->safety;
   stattree->Fill();
   count_gtmedi++;
#endif
}


//______________________________________________________________________
void gmedia(Float_t *x, Int_t &numed, Int_t &check)
{
   fgmedia(x,numed,check);
#ifdef STATISTICS
  statcode = 1;
  statsame = 0;
  for (int j=0;j<6;j++) if (j <3) oldvect[j] = x[j]; else oldvect[j]=0;
  oldsafety = gctrak->safety;
  oldstep   = gctrak->step;
  sprintf(statpath,"%s",geant3->GetPath());
  statsnext=gctrak->snext;
  statsafety=gctrak->safety;
  stattree->Fill();
  count_gmedia++;
#endif
}

//______________________________________________________________________
void gtmany(Int_t &level1)
{
   fgtmany(level1);
}

//______________________________________________________________________
void gtonlyg3(Int_t &isOnly)
{
   //with Geant3, return gonly(nlevel);
   isOnly = (Int_t)gcvolu->gonly[gcvolu->nlevel-1];
}

//______________________________________________________________________
void gtonly(Int_t &isOnly)
{
   //with Geant3, return gonly(nlevel);
   fgtonly(isOnly);
}

//______________________________________________________________________
void glvolu(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier)
{
  //
  //  nlev   number of levels deep into the volume tree
  //         size of the arrays lnam and lnum
  //  lnam   an integer array who's 4 bytes contain the ASCII code for the
  //         volume names
  //  lnum   an integer array containing the copy numbers for that specific
  //         volume
  //
  //  This routine fills the volume parameters in common /gcvolu/ for a
  //  physical tree, specified by the list lnam and lnum of volume names
  //  and numbers, and for all its ascendants up to level 1. This routine
  //  is optimized and does not re-compute the part of the history already
  //  available in GCVOLU. This means that if it is used in user programs
  //  outside the usual framework of the tracking, the user has to initialize
  //  to zero NLEVEL in the common GCVOLU. It return 0 if there were no
  //  problems in make the call.
  //
// printf("glvolu called\n");

  fglvolu(nlev, lnam, lnum, ier);
}


//______________________________________________________________________
void gtnext()
{
#ifdef STATISTICS
   count_gtnext++;
   statcode = 3;
   statsame = 1;
   for (int j=0;j<6;j++) oldvect[j] = gctrak->vect[j];
   oldsafety = gctrak->safety;
   oldstep   = gctrak->step;
   sprintf(statpath,"%s",geant3->GetPath());
#endif

   fgtnext();

#ifdef STATISTICS
  statsnext=gctrak->snext;
  statsafety=gctrak->safety;
  stattree->Fill();
#endif
}
//______________________________________________________________________
void ggperp(Float_t *x, Float_t *norm, Int_t &ierr){
// Computes the normal to the next crossed surface, assuming that
// FindNextBoundary() was already called.

    fggperp(x,norm,ierr);
}
//______________________________________________________________________
Bool_t TGeant3::GetTransformation(const TString &volumePath,TGeoHMatrix &mat){
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
    Int_t i,n,k,*lnam=0,*lnum=0;
    // Default rotation matrix, Unit
    Double_t m[9] = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
    Double_t s[3] = {1.0,1.0,1.0}; // Default scale, Unit
    Double_t t[3] = {0.0,0.0,0.0}; // Default translation, none.

    k =ConvertVolumePathString(volumePath,&lnam,&lnum);//Creates lnam, and lnum
    if(k<=0) { // Error from Convert volumePathString.
        delete[] lnam;
        delete[] lnum;
        return kFALSE;
    } // end if k<=0
    if(k==1){// only one volume listed, must be top most, return unit..
        delete[] lnam;
        delete[] lnum;
        mat.SetRotation(m);
        mat.SetTranslation(t);
        mat.SetScale(s);
        return kTRUE;
    } // end if k==1
    this->Gcvolu()->nlevel = 0;
    i = this->Glvolu(k,lnam,lnum);
    n = this->Gcvolu()->nlevel -1;
    delete[] lnam; // created in ConvertVolumePathString.
    delete[] lnum; // created in ConvertVolumePathString.
    if(i!=0) return kFALSE; // Error
    mat.SetScale(s); // Geant scale always 1.
    if(!((this->Gcvolu()->grmat[n][9])==0.0)) { // not Unit matrix
        for(i=0;i<9;i++) m[i] = (Double_t) this->Gcvolu()->grmat[n][i];
    } // end if
    mat.SetRotation(m);
    for(i=0;i<3;i++) t[i] = (Double_t) (this->Gcvolu()->gtran[n][i]);
    mat.SetTranslation(t);
    return kTRUE;
}/*
//______________________________________________________________________
Bool_t TGeant3::GetShape(const TString &volumeName,TString &shapeType,
                         TArrayD &par){
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
    const Int_t nshapes = 16;
    const Char_t *vname[nshapes] = {"BOX","TRD1","TRD2","TRAP","TUBE","TUBS",
                                   "CONE","CONS","SPHE","PARA","PGON","PCON",
                                   "ELTU","HYPE","GTRA","CTUB"};
    Int_t volid,i,jv0,ishape,npar;
    Float_t *qpar;

    volid = VolId(volumeName.Data());
    if(volid==0) return kFALSE; // Error.
    jv0 = this->Lq()[this->Gclink()->jvolum-volid];
    ishape = (Int_t)(this->Q()[jv0+2]);
    if(ishape<1||ishape>nshapes) return kFALSE; // error unknown shape
    npar   = (Int_t)(this->Q()[jv0+5]);
    qpar   = (this->Q())+jv0+7;
    par.Set(npar); // Resize TArrayD
    for(i=0;i<npar;i++) par.AddAt(((Double_t)qpar[i]),i);
    shapeType = vname[ishape-1];
    return kTRUE;
}*/

//______________________________________________________________________
Bool_t TGeant3::GetShape(const TString &volumePath,TString &shapeType,
                         TArrayD &par){
    // Returns the shape and its parameters for the volume specified
    // by the path volumePath. The format of the path volumePath is as 
    // follows (assuming ALIC is the Top volume)
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
    //   TString &shapeType   Shape type
    //   TArrayD &par         A TArrayD of parameters with all of the
    //                        parameters of the specified shape.
    // Return:
    //   A logical indicating whether there was an error in getting this
    //   information
    const Int_t nshapes = 16;
    const Char_t *vname[nshapes] = {"BOX","TRD1","TRD2","TRAP","TUBE","TUBS",
                                   "CONE","CONS","SPHE","PARA","PGON","PCON",
                                   "ELTU","HYPE","GTRA","CTUB"};
    Int_t volid,i,k,jv0,ishape,npar,*lnam=0,*lnum=0;
    Float_t *qpar;

    k=ConvertVolumePathString(volumePath,&lnam,&lnum);//Creates lnam, and lnum
    if(k<=0) { // Error from Convert volumePathString.
        delete[] lnam;
        delete[] lnum;
        return kFALSE;
    } // end if k<=0
    this->Gcvolu()->nlevel = 0;
    i = this->Glvolu(k,lnam,lnum);
    delete[] lnam;
    delete[] lnum;
    if(i!=0) {// error
        par.Set(0);
        return kFALSE;
    } // end if i!=1
    volid = this->Gcvolu()->lvolum[this->Gcvolu()->nlevel-1];
    jv0 = this->Lq()[this->Gclink()->jvolum-volid];
    ishape = (Int_t)(this->Q()[jv0+2]);
    if(ishape<1||ishape>nshapes) return kFALSE; // error unknown shape
    npar   = (Int_t)(this->Q()[jv0+5]);
    qpar   = (this->Q())+jv0+7;
    par.Set(npar); // Resize TArrayD
    for(i=0;i<npar;i++) par.AddAt((Double_t)qpar[i],i);
    shapeType = vname[ishape-1];
    shapeType = shapeType.Strip();
    return kTRUE;
}
//______________________________________________________________________
Bool_t TGeant3::GetMaterial(const TString &volumeName,
                            TString &name,Int_t &imat,
                            Double_t &a,Double_t &z,Double_t &dens,
                            Double_t &radl,Double_t &inter,TArrayD &par){
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
    Int_t i,volid,jma,nbuf;
    Float_t af,zf,densf,radlf,interf;
    Float_t *ubuf;
    Char_t *ch,namec[20] = {20*'\0'};

    volid = VolId(volumeName.Data());
    if(volid==0) return kFALSE; // Error
    if(volid>0){ // Get Material number, imat.
        Int_t imed = (Int_t) (this->Q()[this->Lq()[
                                            this->Gclink()->jvolum-volid]+4]);
        Int_t jtm  = this->Lq()[this->Gclink()->jtmed-imed];
        imat = (Int_t)(this->Q()[jtm+6]);
    } else {
        i = this->Gclink()->jvolum + volid;
        Int_t jdiv  = this->Lq()[i];
        Int_t ivin  = (Int_t) (this->Q()[jdiv+2]);
        i = this->Gclink()->jvolum - ivin;
        Int_t jvin  = this->Lq()[i];
        Int_t idmed = (Int_t)(this->Q()[jvin+4]);
        i = this->Gclink()->jtmed-idmed;
        Int_t jtm   = this->Lq()[i];
        imat = (Int_t)(this->Q()[jtm+6]);
    } // end if-else
    nbuf = jma = this->Lq()[this->Gclink()->jmate-imat];
    ubuf = new Float_t[nbuf];
    Gfmate(imat,namec,af,zf,densf,radlf,interf,ubuf,nbuf);
    // Problem with getting namec back from Gfmate, get it from 
    // the Zebra bank directly.
    ch = (char *)(this->Iq()+jma+1);
    for(i=0;i<20;i++) if(ch[i]!=' ') namec[i] = ch[i];
    name = namec;
    name = name.Strip();
    //
    par.Set(nbuf);
    for(i=0;i<nbuf;i++) par.AddAt(((Double_t)ubuf[i]),i);
    delete[] ubuf;
    a      = (Double_t) af;
    z      = (Double_t) zf;
    dens   = (Double_t) densf;
    radl   = (Double_t) radlf;
    inter  = (Double_t) interf;
    return kTRUE;
}
//______________________________________________________________________
Bool_t TGeant3::GetMedium(const TString &volumeName,TString &name,
                          Int_t &imed,Int_t &nmat,Int_t &isvol,Int_t &ifield,
                          Double_t &fieldm,Double_t &tmaxfd,Double_t &stemax,
                          Double_t &deemax,Double_t &epsil, Double_t &stmin,
                          TArrayD &par){
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
    Int_t i,volid,nbuf;
    Float_t fieldmf,tmaxfdf,stemaxf,deemaxf,epsilf,stminf;
    Float_t *buf;
    Char_t namec[25] = {25*'\0'};

    volid = VolId(volumeName.Data());
    if(volid==0) return kFALSE; // Error
    if(volid>0){ // Get Material number, imat.
        imed = (Int_t)(this->Q()[this->Lq()[this->Gclink()->jvolum-volid]+4]);
    } else {
        Int_t jdiv  = this->Lq()[this->Gclink()->jvolum + volid];
        Int_t ivin  = (Int_t) (this->Q()[jdiv+2]);
        Int_t jvin  = this->Lq()[this->Gclink()->jvolum - ivin];
        imed = (Int_t)(this->Q()[jvin+4]);
    } // end if-else
    nbuf = this->Lq()[this->Gclink()->jtmed-imed];
    buf  = new Float_t[nbuf];
    Gftmed(imed,namec,nmat,isvol,ifield,fieldmf,tmaxfdf,stemaxf,deemaxf,
           epsilf,stminf,buf,&nbuf);
    name = namec;
    name = name.Strip();
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
//____________________________private method____________________________
Int_t TGeant3::ConvertVolumePathString(const TString &volumePath,
                                       Int_t **lnam,Int_t **lnum){
    // Parses the TString volumePath into an array of volume names 
    // (4 character long limit) and copy numbers in a form used
    // by Geant3.
    // Inputs:
    //   TString& volumePath  The volume path to the specific volume
    //                        for which you want the matrix. Volume name
    //                        hierarchy is separated by "/" while the
    //                        copy number is appended using a "_".
    // Outputs:
    //   Int_t lnam           An integer array, created by this routine,
    //                        containing the 4 character long volume names.
    //   Int_t lnum           An integer array, created by this routine,
    //                        containing the copy numbers.
    // Return:
    //   The size of the arrays lnam an lnum, the number of volumes down
    //   the geometry tree. Note, These arrays are allocated within this
    //   routine, but must be deleted outside of this routine.
    Int_t i,j=0,k=0,ireturn,ichar,*inam;
    Char_t *buf,**levels,**copies,nam[4];

    inam = (Int_t*)nam; // Setup to convert character string to integer.
    buf = new Char_t[volumePath.Length()+1];
    for(i=0;i<volumePath.Length();i++) {
        if(volumePath[i]!=' ')buf[j++] = volumePath[i]; // remove blanks
        if(volumePath[i]=='/') k++;
    } // end for i
    buf[j] = '\0';
    if(buf[j-1]=='/') {k--; buf[j-1]='\0';}// if path ends with '/' ignore 
                                            // it, remove it.
    levels = new Char_t*[k];
    copies = new Char_t*[k];
    (*lnam) = new Int_t[k]; // Allocate Geant3 volume name array
    (*lnum) = new Int_t[k]; // Allocate Geant3 copy number array
    ireturn = k;
    ichar = j;
    k = 0;
    j = 0;
    for(i=0;i<ichar;i++) {
        if(buf[i]=='/'){ 
            levels[k++] = &(buf[i+1]);
            buf[i] = '\0'; // Terminate this sub string.
        } // end if == '/'
        if(buf[i]=='_'){
            copies[j++] = &(buf[i+1]);
            buf[i] = '\0'; // Terminate this sub string.
        } // end if =='_'
    } // end for i
    if(k!=j){ // Error, different number of copy numbers and volume names.
        // clean up everything.
        delete[] buf;
        delete[] levels;
        delete[] copies;
        delete[] (*lnam);
        delete[] (*lnum);
        (*lnam) = 0;
        (*lnum) = 0;
        Error("ConvertVolumePathString","Different number of volume names %d"
              " and copy numbers %d in volumePath:%s",k,j,volumePath.Data());
        return 0;
    } // end if k!=j
    for(i=0;i<k;i++){
        *inam = 0;
        (*lnum)[i] = atoi(copies[i]);
        for(j=0;j<4;j++) {
            if(levels[i][j] == 0) break; // If at end of string exit loop
            nam[j] = levels[i][j];
        } // end for j
        (*lnam)[i] = *inam;
    } // end for i
    // clean up all but lnam and lnum
    delete[] buf;
    delete[] levels;
    delete[] copies;
    return ireturn; // return the size of lnam and lnum. 
}

//__________________________________________________________________
void TGeant3::Gfang( Float_t* p, Float_t& costh, Float_t& sinth, 
                     Float_t& cosph, Float_t& sinph, Int_t& rotate) 
{

  g3fang(p, costh, sinth, cosph, sinph, rotate );
} 

//__________________________________________________________________
Int_t TGeant3::GetIonPdg(Int_t z, Int_t a, Int_t i) const
{
// Acording to
// http://cepa.fnal.gov/psm/stdhep/pdg/montecarlorpp-2006.pdf

  return 1000000000 + 10*1000*z + 10*a + i;
}  
                
//__________________________________________________________________
Int_t TGeant3::GetSpecialPdg(Int_t number) const
{
// Numbering for special particles

  return 50000000 + number;
}                
