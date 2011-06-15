*
* $Id: gdrawc.F 220 2007-11-19 16:08:06Z rdm $
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE G3DRAWC(NAME,IAX,CUTVAL,U0,V0,SU,SV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine draws the intersection of the object        *
C.    *       called NAME with the cut plane normal to one of the      *
C.    *       main axes (IAX=1 is X, IAX=2 is Y, IAX=3 is Z)           *
C.    *       and placed at the distance CUTVAL from the origin.       *
C.    *                                                                *
C.    *       The resulting picture is seen from the axis choosen      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, G3DSPEC                     *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      CHARACTER*4 NAME
      DIMENSION THE(3),PHI(3)
      SAVE THE,PHI
      DATA THE/90.,90.,0./
      DATA PHI/180.,90.,0./
C.
C.    ------------------------------------------------------------------
C.
C             Allocate working space
C
      CALL G3WORK(NWCUT)
C
      ICUT=IAX
      DCUT=CUTVAL
C
      NSURF=0
      ISURF=1
C
C             ICUTFL=1 informs G3DRAW to call GDRWSC instead of G3DRAWS
C             (so that we'll have stored surfaces instead of drawn lines)
C             and to use the line buffer logic of G3DRAWV
C
C             PSI is fixed to 0.
C
      ICUTFL=1
      CALL G3DRAW(NAME,THE(ICUT),PHI(ICUT),0.,U0,V0,SU,SV)
      ICUTFL=0
      END
