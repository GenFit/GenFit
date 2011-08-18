*
* $Id: timel.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: timel.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMVAX)
#include "vaxsys/timel.F"
#else
      SUBROUTINE TIMEL (T)
C
C CERN PROGLIB# Z007    TIMEL   DUMMY   .VERSION KERNFOR  4.11  850529
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      T = 9999.
      RETURN
      END
#endif
