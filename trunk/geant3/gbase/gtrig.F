*
* $Id: gtrig.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gtrig.F,v $
* Revision 1.2  2003/11/28 11:23:55  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3TRIG
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Steering routine to process one event                    *
C.    *    ==>Called by : GRUN   ,<GXINT> GINC4                        *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcvdma.inc"
#include "geant321/gcflag.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
C.
C.    ------------------------------------------------------------------
C.
C.            Kinematics
C.
      IF(RAYTRA.NE.1.)THEN
#if !defined(CERNLIB_USRJMP)
       CALL GUKINE
#endif
#if defined(CERNLIB_USRJMP)
       CALL JUMPT0(JUKINE)
#endif
*      ELSE
*       CALL GDKINE
      ENDIF
      IF (IEOTRI.NE.0)  GO TO  99
C.
C.            Tracking
C.
#if !defined(CERNLIB_USRJMP)
      CALL GUTREV
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUTREV)
#endif
      IF (RAYTRA.EQ.1.) GO TO  99
      IF (IEOTRI.NE.0)  GO TO  99
C.
C.            Digitisation
C.
#if !defined(CERNLIB_USRJMP)
      CALL GUDIGI
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUDIGI)
#endif
      IF (IEOTRI.NE.0)  GO TO  99
C.
C.            User end of event control routine
C.
#if !defined(CERNLIB_USRJMP)
      CALL GUOUT
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUOUT)
#endif
C
  99  RETURN
      END
