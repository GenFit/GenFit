*
* $Id: mzwipe.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mzwipe.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:26  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZWIPE (IXWP)

C-    Wiping, user called

#include "zebra/mqsysh.inc"
#include "zebra/zvfaut.inc"
C--------------    End CDE                             --------------
      DIMENSION    IXWP(9)
#if (defined(CERNLIB_QTRHOLL))&&(!defined(CERNLIB_A6M))
      DIMENSION    NAMESR(2)
      DATA  NAMESR / 4HMZWI, 4HPE   /
#endif
#if (defined(CERNLIB_QTRHOLL))&&(defined(CERNLIB_A6M))
      DATA  NAMESR / 6HMZWIPE /
#endif
#if !defined(CERNLIB_QTRHOLL)
      CHARACTER    NAMESR*8
      PARAMETER   (NAMESR = 'MZWIPE  ')
#endif

#include "zebra/qtrace.inc"
#if defined(CERNLIB_QDEBUG)
      IF (IQVSTA.NE.0)       CALL ZVAUTX
#endif

      IXWIPE = IXWP(1)
      IF (IXWIPE.EQ.0)  IXWIPE=21

      CALL MZGARB (0,IXWIPE)

#include "zebra/qtrace99.inc"
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
