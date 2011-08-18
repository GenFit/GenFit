*
* $Id: qnext.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: qnext.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:47:15  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
#if !defined(CERNLIB_NODUMMY)
      SUBROUTINE QNEXT

C-    DUMMY VERSION OF THE USER ROUTINE QNEXT

#include "zebra/zunit.inc"
C--------------    END CDE                             --------------


      LUN = IQTYPE
      IF (LUN.NE.0)                GO TO 12
   11 LUN = IQLOG
   12 WRITE (LUN,9001)
      IF (LUN.NE.IQLOG)            GO TO 11
      STOP

 9001 FORMAT ('0!!!!!    DUMMY ROUTINE QNEXT REACHED.')
      END
*      ==================================================
#include "zebra/qcardl.inc"
#endif
