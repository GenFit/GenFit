*
* $Id: zabend.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zabend.F,v $
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
      SUBROUTINE ZABEND

C-    DEFAULT FOR THE USER ROUTINE TO HANDLE ABNORMAL PROGRAM END

#include "zebra/zstate.inc"
#include "zebra/zunit.inc"

      CALL ZPOSTM ('TCWM.')
      IF (NQPHAS.GT.0)             GO TO 12
      CALL ABEND
   12 NQPHAS = -2
      CALL ZEND
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
