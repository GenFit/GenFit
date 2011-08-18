*
* $Id: ucopiv.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ucopiv.F,v $
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
      SUBROUTINE UCOPIV (A,B,N)
C
C CERN PROGLIB# V301    UCOPIV          .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71
C
#if defined(CERNLIB_QINTCOPY)
      INTEGER   A(*), B(*)
#endif
#if !defined(CERNLIB_QINTCOPY)
      DIMENSION A(*), B(*)
#endif
C
      IF (N.EQ.0) RETURN
      I = N
      DO 12 J=1,N
      B(I)= A(J)
   12 I = I - 1
      RETURN
      END
