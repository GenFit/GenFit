*
* $Id: gstxin.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gstxin.F,v $
* Revision 1.2  2003/11/28 11:23:56  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE G3STXIN
*   ********************************************************************
*   *                                                                  *
*   *   Initialization of the values needed in computation of          *
*   *   dN^2/dE/dx as in                                               *
*   *   Grishin,Ermilova,Kotelnikov,NIM A307(1991),273                 *
*   *                 PLN 18.04.92, data tables of ILG                 *
*   *                                                                  *
*   *   ==> Called by : G3STINI                                        *
*.  *       Authors   : P. Nevski, K. Lassila-Perini                   *
*   *                                                                  *
*   ********************************************************************
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcstra.inc"
*
      REAL NAV,MB,ME,NE
*
*   1 ev to erg [erg]
      PARAMETER (ERG=1.60218E-12)
*   1mb to cm2
      PARAMETER (MB=1E-18)
*   Electron mass [g]
      PARAMETER (ME=9.10943E-28)
*   Electron charge[ESU]
      PARAMETER (QE=4.80321E-10)
*   Plank constant[erg]
      PARAMETER (H=1.05457E-27)
*
*     ------------------------------------------------------------------
*
*   -----    physic constants
      EPPS=0.01
      EMIN=EM(1)-EPPS
      EMAX=LOG(1.E13)
*     Integration up to infinity
*
      NAV = AVO*1.E24
*
*   Electron mass in ev
      MEEV=ME*CLIGHT**2/ERG
*
*   electron radius[cm]
      R0=QE**2/(ME*CLIGHT**2)
*
*   Plank constant in ev
      HE=H*CLIGHT/ERG
*
*   Gas density
      RO=DENS
*
*   Electron density
      NE=NAV*Z/A*RO
*
*   plasma freq**2 [ev]
      WP2=4*PI*R0*NE*HE**2
*
*   x-section to F.osc
      S1=MB/(2*PI**2*R0*HE*Z)
*
*   dN/dx scale
      S2=2*PI*R0*NE*QE**2/ERG
*
      NTAB = NP
      SFINT=0
*
*    Integration of the oscillator strength function
*    (photoelectric cross-section)
      DO 10  IE=NTAB,1,-1
         E1 = EMM(IE)
         IF (IE.EQ.NTAB) THEN
            E2 = EMAX
         ELSE
            E2 = EMM(IE+1)
         ENDIF
         SFINT = SFINT+ G3OSCIN(EXP(E1),EXP(E2))
         FINT(IE)=SFINT
   10 CONTINUE

*
*   Normalization of the integrated cross-section
      DO 20  IE=1,NTAB
         FINT(IE)=LOG(FINT(IE)/SFINT)
   20 CONTINUE
*
*   Calculation of the real and imaginary parts of the
*   complex dielectric constant
      DO 30  IE=1,NTAB
         E=EXP(EMM(IE))
         EPSI(IE)=LOG(WP2*PIBY2*G3PHSG2(E*1E-9)/(SFINT*E))
*      scaled by
         EPSR(IE)=WP2*G3KOKRI(E,EXP(EMIN),EXP(EMAX))*E/SFINT
   30 CONTINUE
      END
