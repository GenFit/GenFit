*
* $Id: v_313.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: v_313.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:35  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:47  cernlib
* Geant
*
*
#include "sys/CERNLIB_machine.h"
#include "pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.18  by  S.Giani
*-- Author :
*
*::>          VERSION  3.13/05   890628 15.01
*
*        Changes in GTELEC,GTHADR and GTMUON for overstopping tracks.
*        STEP and SLENG correction optimized
*        In COMMON/GCDRAW/ variable IGVIEW changed to IDVIEW
*        Routine GTRMUL recoded
*        DOUBLE precision in GINTCO and protection added.
*        Changes in GNOPG6.
*        Changes in GINME (for spheres,tubes and cones DPP.LE.0)
*        GDINIT called by GXINT
*        Particle data table updated in GPART
*        Tests on boundary conditions for IEKBIN modified in GTGAMA
*        A new user callable routine GBIRK introduced in GPHYS.
*         GBIRK may be called from GUSTEP to compute the Birk
*         correction factors for anorganic scintillators.
*
*::>          VERSION  3.13/04   890623 12.34
*
*        New version of GLANDZ and GPOISS by L.Urban
*        GRNDM calling sequence changed
*        Sequence numbers (1-->215) can be initialized
*         with data card RNDM
*                Ex:
*            RNDM  9876  1234534   initializes sequence 1
*            RNDM 45  initilizes sequence 45 with the starting
*                     seed of sequence 45
*        New routine GRANOR (copy of RANNOR for GRNDM)
*
*::>          VERSION  3.13/03   890425 12.09
*
*        Update of inline documentation DOCGBASE,DOCGKINE,DOCGTRAK
*        Deck GGCLOS modified and moved from GGEOM to GBASE
*        New deck GHCLOS (called by GGCLOS) in GBASE
*        Deck GEVKEV moved from GTRAK to GCONS
*        New decks GFNDIG,GFNHIT,GRHITS in GHITS
*        Minor bug corrected in GTELEC
*        Control of debug and removal of SNXT=BIG  in GTNEXT
*        Call to GUSTEP removed from GUVIEW
*
*::>          VERSION  3.13/02   890311 10.45
*
*        Bug corrected in GPAIRM (COMMON/GCMATE/ overwritten
*        Bug corrected in outines GDRELM and GDRELP
*           (error in AVO changed to AVOGAD)
*        Minor corrections in physics routines
*        SAVE statements added
*        Protection in GHANGL
*        Default changed to GHEISHA instead of TATINA
*          in routines GUPHAD and GUHADR
*        Bug corrected in GHTATI (NAMEC(12) instead of NAMEC(10)
*
