*
* $Id: gtreve.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gtreve.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 07/10/94  18.07.13  by  S.Giani
*-- Author :
      SUBROUTINE G3TREVE
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. G3TREVE                                               *
C.    *                                                                *
C.    *   Controls tracking of all particles belonging to the current  *
C.    *    event.                                                      *
C.    *                                                                *
C.    *   Called by : GUTREV, called by G3TRIG                         *
C.    *   Authors   : R.Brun, F.Bruyant                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
      REAL UBUF(2)
      EQUIVALENCE (UBUF(1),WS(1))
      LOGICAL   BTEST
C.
C.    ------------------------------------------------------------------
      NTMSTO = 0
      NSTMAX = 0
      NALIVE = 0
*
      NV = NVERTX
      DO 290 IV = 1,NV
*
* ***   For each vertex in turn ..
*
         JV = LQ(JVERTX-IV)
         NT = Q(JV+7)
         IF (NT.LE.0) GO TO 290
         TOFG   = Q(JV+4)
         SAFETY = 0.
*
*  **   Loop over tracks attached to current vertex
*
         DO 190 IT = 1,NT
            JV   = LQ(JVERTX-IV)
            ITRA = Q(JV+7+IT)
            IF (BTEST(IQ(LQ(JKINE-ITRA)),0)) GO TO 190
            CALL G3FKINE (ITRA, VERT, PVERT, IPART, IVERT, UBUF, NWBUF)
            IF (IVERT.NE.IV) THEN
               WRITE (CHMAIL, 1002) IV, IVERT
               CALL GMAIL (0, 0)
               GO TO 999
            ENDIF
*
*   *      Store current track parameters in stack JSTAK
*
            CALL G3SSTAK (2)
  190    CONTINUE
*
*  **   Start tracking phase
*
  210    IF (NALIVE.NE.0) THEN
            NALIVE = NALIVE -1
*
*   *      Pick-up next track in stack JSTAK, if any
*
            IF (IQ(JSTAK+1).GT.0) THEN
*
*   *         Initialize tracking parameters
*
               CALL G3LTRAC
               IF (NUMED.EQ.0) GO TO 210
            ELSE
*
*   *         otherwise, select next track segment from stack JTRACK
*
               print *,' =======CALLING GFTRAC'
c               CALL GFTRAC
*
            ENDIF
*
*   *       Resume tracking
*
#if !defined(CERNLIB_USRJMP)
              CALL GUTRAK
#endif
#if defined(CERNLIB_USRJMP)
              CALL JUMPT0(JUTRAK)
#endif
            IF (IEOTRI.NE.0) GO TO 999
            GO TO 210
         ENDIF
*
  290 CONTINUE
*
 1001 FORMAT (' GTREVE : Vertex outside setup, XYZ=',3G12.4)
 1002 FORMAT (' GTREVE : Abnormal track/vertex connection',2I8)
*                                                             END GTREVE
  999 END
 
