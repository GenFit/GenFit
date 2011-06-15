*
* $Id: gdrawx.F 220 2007-11-19 16:08:06Z rdm $
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE G3DRAWX(NAME,CUTTHE,CUTPHI,CUTVAL,THE,PHI,U0,V0,SU,SV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw extended cut : this routine is similar to G3DRAWC   *
C.    *       but any cut plane can be choosen, instead of just        *
C.    *       the three planes normal to the main axes.                *
C.    *                                                                *
C.    *       The cut plane is normal to the line given by             *
C.    *       the two cut angles CUTTHE,CUTPHI                         *
C.    *       and it is placed at the distance CUTVAL from the origin. *
C.    *                                                                *
C.    *       The resulting picture is seen from the angles THE,PHI    *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      CHARACTER*4 NAME
C.
C.    ------------------------------------------------------------------
C.
C             Allocate working space
C
      CALL G3WORK(NWCUT)
C
      IF (CUTTHE.EQ.90.AND.CUTPHI.EQ.180) THEN
         ICUT=1
      ELSE IF (CUTTHE.EQ.90.AND.CUTPHI.EQ.0) THEN
         ICUT=1
      ELSE IF (CUTTHE.EQ.90.AND.CUTPHI.EQ.90) THEN
         ICUT=2
      ELSE IF (CUTTHE.EQ.0.AND.CUTPHI.EQ.0) THEN
         ICUT=3
      ELSE
         ICUT=0
      ENDIF
C
      DCUT=CUTVAL
      CTHETA=CUTTHE
      CPHI=CUTPHI
C
      NSURF=0
      ISURF=1
C
C             ICUTFL=1 informs G3DRAW to call GDRWSC instead of G3DRAWS
C             (so that we'll have stored surfaces instead of drawn lines
C             and to use the line buffer logic of G3DRAWV
C
C             PSI is fixed to 0.
C
      ICUTFL=1
      CALL G3DRAW(NAME,THE,PHI,0.,U0,V0,SU,SV)
      ICUTFL=0
      END
