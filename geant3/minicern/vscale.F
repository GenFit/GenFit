*
* $Id: vscale.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: vscale.F,v $
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
      SUBROUTINE VSCALE (B,ALPHA,A,N)
C
C CERN PROGLIB# F121    VSCALE          .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION A(*),B(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= B(I) * ALPHA
C
      RETURN
      END
