*
* $Id: iufind.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: iufind.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:11  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION IUFIND (IT,ICH,JL,JR)
C
C CERN PROGLIB# V304    IUFIND          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION ICH(9)
C
      DO 12 J= JL,JR
      IF (ICH(J).EQ.IT)   GO TO 24
   12 CONTINUE
      J= JR+1
C
   24 IUFIND= J
      RETURN
      END
