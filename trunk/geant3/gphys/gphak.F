*
* $Id: gphak.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gphak.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      REAL FUNCTION G3PHAK(BET)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates Photoelectron Angular ditribution from K shell      *
C.    *   from L1 shell ditribution ENTRY G3PHAL1                      *
C.    *                                                                *
C.    *    ==>CALLED BY : G3PHOT                                       *
C.    *       AUTHOR    : J. Chwastowski                               *
C.    *                                                                *
C.    ******************************************************************
C.
      DOUBLE PRECISION SAUT1,SAUT2,SAUT3,DRN
      DOUBLE PRECISION COST,SIN2T,X,ANOR,AMAJOR
      DOUBLE PRECISION GAMA,G1,G2,BETA,B1,ONE,HALF
      PARAMETER (BETMAX=0.999999)
      PARAMETER (ONE=1,HALF=ONE/2)
      DIMENSION RAN(2),BETVAL(7),AMAJOR(7)
      DATA BETVAL / 0.8601,0.901,0.9901,0.99901,0.999901,0.9999901,
     +              0.99999901 /
      DATA AMAJOR / 1.02,1.1,4.3,19.,68.,221.,704. /
      DATA ANOR / 0.0 /
C
C L1 shell entry
C
      ENTRY G3PHAL1(BET)
C
      BETA = MIN(BET,BETMAX)
      B1 = ONE-BETA
      G2 = ONE/((ONE-BETA)*(ONE+BETA))
      GAMA = SQRT(G2)
      G1 = GAMA-ONE
      DO 10 I = 1,7
       IF(BETA.LE.BETVAL(I)) THEN
        ANOR = AMAJOR(I)*G2
        ANOR = ONE/ANOR
        GO TO 20
       ENDIF
   10 CONTINUE
   20 SAUT2 = HALF*GAMA*G1*(GAMA-2)
      G3PHAK = 0.0
   30 CONTINUE
      CALL GRNDM(RAN,2)
      DRN = 2*RAN(1)
      SAUT1 = (DRN-B1)/(DRN*BETA+B1)
      COST = SAUT1
      X = ONE-BETA*COST
      SIN2T = (ONE-COST)*(ONE+COST)
      SAUT3 = SIN2T/(X**2)*(ONE+SAUT2*X)*ANOR
      IF(RAN(2).GT.SAUT3) GO TO 30
      G3PHAK  = COST
      G3PHAL1 = G3PHAK
      END
