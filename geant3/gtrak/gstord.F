*
* $Id: gstord.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gstord.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE G3STORD
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Orders the particles in the stack in ascending order by      *
C.    *   the user word.                                               *
C.    *                                                                *
C.    *   Called by : G3LTRAC                                          *
C.    *   Authors   : F.Carminati                                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcnum.inc"
      PARAMETER (NWPOLA=3,ISTPOS=12)
*
* ----------------------------------------------------------------------
*
      JSTPOS = JSTAK+3+ISTPOS
      VALMAX = Q(JSTPOS)
      INDMAX = 1
      NPSTAK = IQ(JSTAK+1)
      DO 10 J=2,IQ(JSTAK+1)
         JSTPOS = JSTPOS+NWSTAK
         VALTST = Q(JSTPOS)
         IF(VALTST.GT.VALMAX) THEN
            VALMAX = VALTST
            INDMAX = J
         ENDIF
   10 CONTINUE
      IF(INDMAX.NE.NPSTAK) THEN
         JSMOVE = JSTAK+3+(INDMAX-1)*NWSTAK
         JSLAST = JSTAK+3+(NPSTAK-1)*NWSTAK
         DO 20 J=1,NWSTAK
            TEMP = Q(JSMOVE+J)
            Q(JSMOVE+J) = Q(JSLAST+J)
            Q(JSLAST+J) = TEMP
   20    CONTINUE
         IF(JPOLAR.NE.0) THEN
            JPMOVE = JPOLAR+(INDMAX-1)*NWPOLA
            JPLAST = JPOLAR+(NPSTAK-1)*NWPOLA
            DO 30 J=1,NWPOLA
               TEMP = Q(JPMOVE+J)
               Q(JPMOVE+J) = Q(JPLAST+J)
               Q(JPLAST+J) = TEMP
   30       CONTINUE
        ENDIF
      ENDIF
      END
