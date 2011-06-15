*
* $Id: gfstat.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfstat.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFSTAT(ICASE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       FILL banks for volume statistics                         *
C.    *                                                                *
C.    *    ==>Called by : GINVOL, GMEDIA, GTNEXT, GTMEDI, GTRACK       *
C.    *                                                                *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF(JGSTAT.LE.0)GO TO 999
      IF(NLEVEL.LE.0)GO TO 999
C
      IF(ICASE.LE.1) THEN
C
C                 GINVOL   ICASE=0,1
C
         IVO=LVOLUM(NLEVEL)
         JG=LQ(JGSTAT-IVO)
         IF(JG.GT.0) THEN
            IQ(JG+1+ICASE)=IQ(JG+1+ICASE)+1
         ENDIF
      ELSEIF(ICASE.LE.10) THEN
C
C         GMEDIA  (ICASE=2), GTNEXT (ICASE=3) or GTMEDI (ICASE=4)
C
         IC = ICASE*2-1
         DO 10 ILEV=1,NLEVEL
            IVO=LVOLUM(ILEV)
            JG=LQ(JGSTAT-IVO)
            IF(JG.GT.0) THEN
               IQ(JG+IC)=IQ(JG+IC)+1
               IF(ILEV.EQ.NLEVEL) THEN
                  IQ(JG+IC+1)=IQ(JG+IC+1)+1
               ENDIF
            ENDIF
   10    CONTINUE
C
C               ICASE=10+ITRTYP  (steps)
C
      ELSE
         IVO=LVOLUM(NLEVEL)
         JG=LQ(JGSTAT-IVO)
         IF(JG.GT.0) THEN
            IQ(JG+9)=IQ(JG+9)+1
         ENDIF
      ENDIF
C
  999 END
