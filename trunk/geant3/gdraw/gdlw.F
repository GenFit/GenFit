*
* $Id: gdlw.F 220 2007-11-19 16:08:06Z rdm $
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDLW(LW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Set line width to LW (between 1 to 5) :                  *
C.    *       if LW > 0 set line width temporarily                     *
C.    *       if LW < 0 set line width permanently                     *
C.    *       if LW = 0 restore permanent line width value             *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, G3DHEAD, G3DTREE            *
C.    *       Author : P.Zanarini ; S.Giani  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
***      LINW=MOD(IABS(LW),8)
      IF (LW.EQ.0) THEN
         CALL MVBITS(LINATP,7,3,LINATT,7)
      ELSE IF (LW.GT.0) THEN
         CALL MVBITS(LW,0,3,LINATT,7)
      ELSE
         CALL MVBITS(LW,0,3,LINATT,7)
         CALL MVBITS(LW,0,3,LINATP,7)
      ENDIF
      END
