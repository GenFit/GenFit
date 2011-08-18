*
* $Id: gtmed.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gtmed.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:43  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE GTMED(MEDGEA,MEDMOR)
#include "geant321/mmicap.inc"
#include "geant321/minput.inc"
C get MICAP material number
      DO 10 I=1,MEDIA
         IF(LD(LGE2MO+I).EQ.MEDGEA) THEN
            MEDMOR = I
            GOTO 20
         ENDIF
   10 CONTINUE
      WRITE(IOUT,'('' MICAP GTMED: GEANT Medium '',I5,
     +             '' not found ==> STOP'')') MEDGEA
      STOP
   20 RETURN
      END
