*
* $Id: gbteth.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gbteth.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      REAL FUNCTION G3BTETH(ENER,PARTM,EFRAC)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Universal distribution function for the Bremsstrahlung       *
C.    *   and pair production angles. This function approximates       *
C.    *   the real distribution function which can be found in         *
C.    *   Yung-Su Tsai: Rev. Mod. Phys. 46,815(1974)                   *
C.    *     +erratum  : Rev. Mod. Phys. 49,421(1977)                   *
C.    *                                                                *
C.    *   Called by : G3BREME, G3BREMM, G3PAIRG                        *
C.    *   Authors   : L.Urban                                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
      DIMENSION RNDM(3)
      PARAMETER (ALFA=0.625)
*     ----------------------------------------------------------
* GEANT 3.15
*     G3BTETH=1.
*     GOTO 999
* Angular distribution with d=const
*     D=27.
* Angular distribution with d=d(Z,E,y)
      D=0.13*(0.8+1.3/Z)*(100.+1./ENER)*(1.+EFRAC)
*     -------------------------------------------------------------
*
      W1=9./(9.+D)
      UMAX=ENER*PI/PARTM
*
10    CALL GRNDM(RNDM,3)
*
      IF(RNDM(1).LE.W1) THEN
         BETA=ALFA
      ELSE
         BETA=3.*ALFA
      ENDIF
*
      U=-(LOG(RNDM(2)*RNDM(3)))/BETA
*
*     cut : theta should be .LE. PI  !!
*     this condition depend on E in the case of d=const too!!!!!
*
      IF(U.GE.UMAX) GOTO 10
*
      G3BTETH=U
*
 999  END
