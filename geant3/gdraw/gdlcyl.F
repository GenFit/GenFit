*
* $Id: gdlcyl.F 220 2007-11-19 16:08:06Z rdm $
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDLCYL(R1,PHI1,Z1,R2,PHI2,Z2)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw a line between points (R1,PHI1,Z1) and (R2,PHI2,Z2) *
C.    *       given in cylindrical coordinates                         *
C.    *                                                                *
C.    *    ==>Called by : G3DRAWS                                      *
C.    *       Author : A.McPherson   *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
      DIMENSION P1(3),P2(3)
C.
C.    ------------------------------------------------------------------
C
      IF (R1.LT.0.0001.AND.R2.LT.0.0001) GO TO 999
C
      PH=PHI1*DEGRAD
      P1(1)=R1*COS(PH)
      P1(2)=R1*SIN(PH)
      P1(3)=Z1
      PH=PHI2*DEGRAD
      P2(1)=R2*COS(PH)
      P2(2)=R2*SIN(PH)
      P2(3)=Z2
C
      CALL GDLINE(P1,P2)
C
  999 RETURN
      END
