*
* $Id: zend.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zend.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:26  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:47:14  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE ZEND

C-    DEFAULT FOR THE USER ROUTINE FOR PROGRAM STOP

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"
C--------------    END CDE                             --------------


      CALL ZPHASE (-1)
      CALL MZEND

      IF (NQERR.EQ.0)              GO TO 99
      CALL ABEND
   99 CONTINUE
      STOP
      END
*      ==================================================
#include "zebra/qcardl.inc"
