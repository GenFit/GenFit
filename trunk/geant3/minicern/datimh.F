*
* $Id: datimh.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: datimh.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
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
#if defined(CERNLIB_QMMPW)
#include "mpwgs/datimh.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/datimh.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/datimh.F"
#else
      SUBROUTINE DATIMH (ND,NT)
C
C CERN PROGLIB# Z007    DATIMH  DUMMY   .VERSION KERNFOR  4.03  821008
C
C-    DUMMY FOR NON-ESSENTIAL ROUTINE STILL MISSING ON YOUR MACHINE

      DIMENSION ND(9), NT(9)
      DIMENSION M(8)

#if defined(__GNUC__) && __GNUC__ >= 4
      CALL UBLOW  ('29/09/79', M, 8)
      CALL UBUNCH (M, ND, 8)
      CALL UBLOW  ('12.00.00',M,8)
      CALL UBUNCH (M,NT,8)
#else
      CALL UBLOW (8H29/09/79,M,8)
      CALL UBUNCH           (M,ND,8)
      CALL UBLOW (8H12.00.00,M,8)
      CALL UBUNCH           (M,NT,8)
#endif
      RETURN
      END
#endif
