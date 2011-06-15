*
* $Id: iusame.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: iusame.F,v $
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
      FUNCTION IUSAME (IA,JST,N,MIN,JSAME)
C
C CERN PROGLIB# M501    IUSAME          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71 CL
C
      DIMENSION IA(*)
C
C
      J  = JST
   21 NSM= 0
      JSM= J
C
   24 NSM= NSM + 1
      IF (J.EQ.N)  GO TO 31
      J = J + 1
      IF (IA(J).EQ.IA(J-1))  GO TO 24
C
      IF (NSM.LT.MIN)  GO TO 21
      GO TO 32
C
   31 IF (NSM.GE.MIN)  GO TO 32
      NSM= 0
      JSM= N + 1
   32 IUSAME=NSM
      JSAME= JSM
      RETURN
      END
