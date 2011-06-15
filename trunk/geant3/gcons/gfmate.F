*
* $Id: gfmate.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfmate.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3FMATE(IMAT,NAMATE,A,Z,DENS,RADL,ABSL,UBUF,NWBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Return parameters for material IMAT                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GTRACK                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER NAMATE*(*)
      DIMENSION UBUF(1)
C.
C.    ------------------------------------------------------------------
C.
      A=-1.
      IF (JMATE.LE.0)GO TO 99
      IF (IMAT.LE.0)GO TO 99
      IF (IMAT.GT.NMATE)GO TO 99
C
      JMA = LQ(JMATE- IMAT)
      IF (JMA.LE.0)GO TO 99
      CALL UHTOC(IQ(JMA+1),4,NAMATE,MIN(20,LEN(NAMATE)))
      A = Q(JMA + 6)
      Z = Q(JMA + 7)
      DENS = Q(JMA + 8)
      RADL = Q(JMA + 9)
      ABSL = Q(JMA + 10)
      NWBUF = IQ(JMA-1) - 11
      IF(NWBUF.GT.0) CALL UCOPY(Q(JMA+12),UBUF,NWBUF)
C
  99  RETURN
      END
