*
* $Id: gfrung.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfrung.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3FRUNG (NWRUNG,IRUNG,NWBUF,UBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C     *       Retrieves the parameter of the RUNG   bank               *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *       Author     M.Maire                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
*
      DIMENSION IRUNG(1),UBUF(1)
*
      IF (JRUNG.NE.0) THEN
         NWRUNG = IQ(JRUNG-1)
         DO 10 I=1,NWRUNG
            IRUNG(I) = IQ(JRUNG+I)
   10    CONTINUE
*
         NWBUF = 0
         IF (LQ(JRUNG-1).NE.0) THEN
            JRU = LQ(JRUNG-1)
            NWBUF = IQ(JRU-1)
            DO 20 I=1,NWBUF
               UBUF(I) = Q(JRU+I)
   20       CONTINUE
         ENDIF
*
      ELSE
         NWRUNG = 0
         WRITE(CHMAIL,1000)
         CALL GMAIL(0,0)
      ENDIF
*
 1000 FORMAT(' **** GFRUNG: No RUNG bank existing in memory ')
*
      END
