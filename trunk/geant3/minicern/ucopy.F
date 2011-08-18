*
* $Id: ucopy.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ucopy.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:19  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UCOPY (A,B,N)
C
C CERN PROGLIB# V301    UCOPY           .VERSION KERNFOR  4.40  940929
C ORIG. 01/01/65 JZ
C
#if defined(CERNLIB_QINTCOPY)
      INTEGER   A(*),B(*)
#endif
#if !defined(CERNLIB_QINTCOPY)
      DIMENSION A(*),B(*)
#endif
C
C--                NO OVERLAP  OR  BEGINNING OF A ON END OF B
      IF (N.EQ.0) RETURN
         DO 21 I=1,N
   21 B(I)=A(I)
      RETURN
      END
