*
* $Id: ufill.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ufill.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UFILL (IV,JL,JR,ISTUFF)
C
C CERN PROGLIB# V300    UFILL           .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71 JZ
C
      DIMENSION  IV(*)
C
      ITEXT=ISTUFF
         DO 29 J=JL,JR
   29 IV(J)=ITEXT
      RETURN
      END
