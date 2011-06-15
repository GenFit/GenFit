*
* $Id: gbstat.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gbstat.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE G3BSTAT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Book banks for volume statistics                         *
C.    *                                                                *
C.    *    ==>Called by : G3GCLOS                                      *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF(JGSTAT.GT.0)GO TO 999
      IF(JVOLUM.LE.0)GO TO 999
      CALL MZBOOK(IXCONS,JGSTAT,JGSTAT,1,'GSTA',NVOLUM,NVOLUM,0,2,0)
*
      DO 10 I=1,NVOLUM
         CALL MZBOOK(IXCONS,LBANK,JGSTAT,-I,'GSTA',0,0,9,2,0)
   10 CONTINUE
*
  999 END
