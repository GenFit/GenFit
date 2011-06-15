*
* $Id: gfnrdt.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfnrdt.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:25  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE G3FNRDT(Z,ESHL,NSHL,NWORD,NONRAD,PRB,ENER)
C.
C.    ******************************************************************
C.    *                                                                *
C     *   Fetch Non-Radiative Decay probability Table and the modes    *
C.    *       the nonradiative shells decay                            *
C.    *                                                                *
C.    *    ==>Called by : G3FSHDC                                      *
C.    *       Author    J. Chwastowski                                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gccuts.inc"
#include "geant321/gcpmxz.inc"
#include "geant321/gcphnr.inc"
      REAL  ONEEV
      PARAMETER (ONEEV = 1.E-9)
      DIMENSION ESHL(*),PRB(*),ENER(*),NONRAD(*)
      IZ = Z
      NWORD = 0
      DO 10 J = 1,NSHL
         NONRAD(J) = 0
         ESHELL = ESHL(J)
         IF(ESHELL.GT.0.0) THEN
            I1 = IGNRFN((J-1)*2+1,IZ)
            I2 = IGNRFN((J-1)*2+2,IZ)
            IF(I1*I2.GT.0) THEN
               NONRAD(J) = 1
               NWORD = NWORD+1
               ENER(NWORD) = ESHELL-ESHL(I1)-ESHL(I2)
               ENER(NWORD) = ENER(NWORD)*ONEEV
               IF(ENER(NWORD).LE.CUTELE) ENER(NWORD) = -ENER(NWORD)
               PRB(NWORD) = 1.
            ENDIF
         ENDIF
   10 CONTINUE
      END
