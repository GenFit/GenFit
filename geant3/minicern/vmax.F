*
* $Id: vmax.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: vmax.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:17  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION VMAX (A,N)
C
C CERN PROGLIB# F121    VMAX            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*)
C
      IF (N.LE.0)  GO TO 100
      FMAX = A(1)
      DO 9 I= 1,N
      FMAX = MAX (A(I),FMAX)
    9 CONTINUE
      VMAX = FMAX
      RETURN
  100 VMAX = 0.
      RETURN
      END
