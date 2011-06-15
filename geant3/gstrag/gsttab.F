*
* $Id: gsttab.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gsttab.F,v $
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
      SUBROUTINE G3STTAB(GAMMA,NT,EN,FN)
C   ********************************************************************
C   *                                                                  *
C   *  dN^2/dEdx table for a gamma factor                              *
C   *                                                                  *
C   *  ==> Called by : G3STINI                                         *
C   *      Authors   : P. Nevski                                       *
C   *                                                                  *
C   ********************************************************************
#include "geant321/gcstra.inc"
      EXTERNAL G3STDN
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION DGAMM, DGAM2, DBET2
#endif
 
      DIMENSION EN(1),FN(1)
      PARAMETER (RENER=0.05)
*
*     ------------------------------------------------------------------
*
C   ----  final tables for dN/dx integral vs Energy for a gamma   ------
      DGAMM = GAMMA
      DGAM2 = DGAMM**2
      GAMMA2= DGAM2
      DBET2 = 1-1/DGAM2
      BETA2 = DBET2
      S=0
      NT=0
      ES=0
      EI=EMAX
C   Loop over I=NTAB,1,-1
      DO 10 I=NTAB,1,-1
         EO=EI
         EI=EMM(I)
         S=S+G3XGINT(G3STDN,EI,EO,EPPS)
C      Skip Unless ABS(ES-EI)>RENER \ I==1
         IF ((ABS(ES-EI).LE.RENER).AND.(I.NE.1))GO TO 10
         ES=EI
         NT=NT+1
         EN(NT)=EI
         FN(NT)=S
   10 CONTINUE
      DNDX=S
      END
