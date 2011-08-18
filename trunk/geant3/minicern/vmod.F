*
* $Id: vmod.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: vmod.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:18  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VMOD (A,N)
C
C CERN PROGLIB# F121    VMOD            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      XX= 0.
      IF (N.LE.0)  GO TO 100
      DO 9 I= 1,N
    9 XX= XX + A(I)*A(I)
      IF (XX.LE.0.)  GO TO 100
      VMOD= SQRT (XX)
      RETURN
  100 VMOD = 0.
      RETURN
      END
