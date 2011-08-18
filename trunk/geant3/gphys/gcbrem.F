*
* $Id: gcbrem.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcbrem.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      FUNCTION G3CBREM(T)
C.    ******************************************************************
C.    *                                                                *
C.    *  It calculates the ratio of positron to electron               *
C.    *       bremsstrahlung energy loss.                              *
C.    *                                                                *
C.    *  Input:kinetic energy in GeV                                   *
C.    *  Output: ratio of positron to electron brems energy loss       *
C.    *                                                                *
C.    *  ( The ratio is calculated from a function fitted to the       *
C.    *    results of L.Kim et al.,PHYS.REV.A33,3002(1986)    )        *
C.    *                                                                *
C.    *  Author   L. Urban                                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
      PARAMETER (AA=7522100.,A1=0.415,A3=0.0021,A5=0.00054)
C.
      CFAC=0.
      IF(Z.LE.0.) GOTO 10
      X=LOG(AA*T/Z**2)
*
      IF(X.LE.-8.) GOTO 10
*
      IF(X.GE.+9.) THEN
         CFAC=1.
      ELSE
         W=A1*X+A3*X**3+A5*X**5
         CFAC=0.5+ATAN(W)/PI
      ENDIF
*
   10 G3CBREM = CFAC
*
      END
