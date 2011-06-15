*
* $Id: uoptc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: uoptc.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:51  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UOPTC (DATA,POSS,IOPTV)
C
C CERN PROGLIB# M502    UOPTC           .VERSION KERNFOR  4.10  850320
C ORIG. 07/11/84, JZ
C
      INTEGER      IOPTV(9)
      CHARACTER    *(*)  DATA, POSS


      NP = LEN(POSS)

      DO 24 J=1,NP
      IOPTV(J) = 0
      IF (INDEX(DATA,POSS(J:J)).NE.0)  IOPTV(J)= 1
   24 CONTINUE
      RETURN
      END
