*
* $Id: glast.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: glast.F,v $
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
      SUBROUTINE G3LAST
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       TERMINATION ROUTINE                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGLAST                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctime.inc"
C.
C.    ------------------------------------------------------------------
C.
      WRITE (CHMAIL,1000) IEVENT
      CALL GMAIL(0,0)
C
      CALL GRNDMQ(NRNDM(1),NRNDM(2),0,'G')
C
      WRITE (CHMAIL,3000) NRNDM
      CALL GMAIL(0,0)
C
C             COMPUTE ONE EVENT PROCESSING TIME
C
      IF(IEVENT.GT.0)THEN
         CALL TIMEL(TIMLFT)
         XMEAN  = (TIMINT - TIMLFT)/IEVENT
         WRITE(CHMAIL,4000)XMEAN
         CALL GMAIL(0,2)
      ENDIF
C
C             Print ZEBRA statistics
C
      CALL MZEND
C
C             Print statistics
C
      IF(JGSTAT.NE.0)CALL G3PSTAT
#if defined(CERNLIB_MONITOR)
C
C             Call monitoring routine
C
      CALL G3EAMON(2,'GLAST')
#endif
C
 1000 FORMAT('1',9X,'**** NUMBER OF EVENTS PROCESSED =',I10)
 3000 FORMAT(10X,'**** RANDOM NUMBER GENERATOR AFTER'
     +,' LAST COMPLETE EVENT ',2I12)
 4000 FORMAT(10X,'**** TIME TO PROCESS ONE EVENT IS =',F10.4,
     +      ' SECONDS')
      END
