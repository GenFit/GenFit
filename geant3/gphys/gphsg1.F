*
* $Id: gphsg1.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gphsg1.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
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
      FUNCTION G3PHSG1(E)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       FUNCTION TO COMPUTE THE PHOTOEFFECT TOTAL CROSS-SECTION  *
C.    *       AS A FUNCTION OF E                                       *
C.    *                                                                *
C.    *    ==>CALLED BY : G3PHOTI                                      *
C.    *                                                                *
C.    *    (1) Sandia parametrizsation is used;                        *
C.    *    (2) GPHSG1 in cm**2/g.                                      *
C.    *    (3) E in GeV                                                *
C.    *                                                                *
C.    *     J. Chwastowski 17.11.92                                    *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gc10ev.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION EINV,ECUR,ONE
#endif
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      RES = 0.
      IF(E.LT.G10EV) GO TO 20
C Use Sandia data
      JPHXS = LQ(JPHOT-1)
      IPOINT = JPHXS+Q(JPHXS+1)*3+2
      IMAX = Q(IPOINT)
      IPOINT = IPOINT+1
      ECUR = E*1.E6
      IF(ECUR.LT.Q(IPOINT)) GO TO 20
      EINV = ONE/ECUR
      DO 10 I = 2,IMAX
         IPOINT = IPOINT+5
         IF(ECUR.LT.Q(IPOINT)) THEN
            J = IPOINT+1
            RES = EINV*(Q(J)+EINV*(Q(J+1)+EINV*(Q(J+2)+EINV*Q(J+3))))
            GO TO 20
         ENDIF
   10 CONTINUE
C This value is in cm**2/g
   20 G3PHSG1 = MAX(RES,0.)
C
      END
