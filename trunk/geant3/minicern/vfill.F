*
* $Id: vfill.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: vfill.F,v $
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
      SUBROUTINE VFILL (A,N,STUFF)
C
C CERN PROGLIB# F121    VFILL           .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71
C
#if defined(CERNLIB_QINTCOPY)
      INTEGER   A(*), STUFF
#endif
#if !defined(CERNLIB_QINTCOPY)
      DIMENSION A(*)
#endif
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= STUFF
C
      RETURN
      END
