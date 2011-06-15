*
* $Id: gmiset.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gmiset.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:43  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  18.26.26  by  S.Giani
*-- Author :    S.Ravndal   22/02/95
       SUBROUTINE GMISET(MATNO,NKEY)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       MIFACE Micap Interface routine to call MICSET            *
C.    *                                                                *
C.    *   ==> Called by : <User>                                       *
C.    *       Authors   : S.Ravndal                                    *
C.    *                                                                *
C     *       for details see MICAP manual ORNL/TM-10340               *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
 
       INTEGER MATNO,NKEY
C
       IF (NKEY.GE.0.AND.NKEY.LE.1) THEN
              CALL MICSET(MATNO,NKEY,' ')
       ELSE
              WRITE(CHMAIL,10100)
              CALL GMAIL(0,0)
       END IF
       RETURN
10100  FORMAT(1X,'GMISET called with non-valid parameters')
       END
 
