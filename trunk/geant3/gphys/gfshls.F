*
* $Id: gfshls.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfshls.F,v $
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
      SUBROUTINE G3FSHLS(Z,GSHPOT,NSHELL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Fetch the shells' potentials for e given Z               *
C.    *       The potentials are in eV !!!                             *
C.    *                                                                *
C.    *    ==>Called by : G3PHXSI,G3FSHDC                              *
C.    *       Author    J. Chwastowski                                 *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT NONE
#include "geant321/gcpmxz.inc"
#include "geant321/gcshpt.inc"
#include "geant321/gcunit.inc"
      REAL ZZ,Z,GSHPOT(*)
      INTEGER I,IZ,NSHELL
      ZZ = Z
      IF(ZZ.LT.1.OR.ZZ.GT.MAXELZ) THEN
         IF(ZZ.LT. 1) WRITE(CHMAIL,10000) ZZ
         IF(ZZ.GT.MAXELZ) WRITE(CHMAIL,10100) ZZ
         CALL GMAIL(0,0)
      ELSE
         IZ = IFIX(ZZ)
         NSHELL = NSHLLS(IZ)
         DO 10 I = 1,NSHLLS(IZ)
            GSHPOT(I) = ESHELL(N1ST(IZ)-1+I)
   10    CONTINUE
      ENDIF
10000 FORMAT(' ***** GFSHLS ERROR:',
     +       '  Z of the material is less than   1.',
     +       ' Actual Z =',I5,'.')
10100 FORMAT(' ***** GFSHLS ERROR:',
     +       '  Z of the material is bigger than 100.',
     +       ' Actual Z =',I5,' (I5).')
      END
