*
* $Id: traceq.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: traceq.F,v $
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
#if defined(CERNLIB_QMHPX)
#include "hpxgs/traceq.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxsys/traceq.F"
#else
      SUBROUTINE TRACEQ (LUNP,LEV)
C
C CERN PROGLIB# N105    TRACEQ          .VERSION KERNFOR  4.23  891215
C ORIG.  4/01/89  JZ
C
*     unit number for standard output
      PARAMETER    (LUNSTP=6)

      LUN = LUNP
      IF (LUN.EQ.0)  LUN = LUNSTP
      WRITE (LUN,9001)
      RETURN
 9001 FORMAT (/' TRACEQ.  In-line trace-back still not available.')
      END
#endif
