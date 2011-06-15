*
* $Id: mzgsta.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mzgsta.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:23  fca
* AliRoot sources
*
* Revision 1.2  1996/04/18 16:12:36  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:47:20  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZGSTA (IGARB)

C-    Bump garbage collection statistic

#include "zebra/mqsys.inc"
#include "zebra/mzct.inc"
C--------------    End CDE                             --------------
      DIMENSION    IGARB(20)


      LMT  = LQMTA
   22 IACT = LQ(LMT+1)
      IF (IACT.LT.3)               GO TO 28
      JDIV = LQ(LMT)
      IF (IACT.EQ.3)               GO TO 26
      NQDWIP(KQT+JDIV) = NQDWIP(KQT+JDIV) + 1
      GO TO 28

   26 IGARB(JDIV) = IGARB(JDIV) + 1
   28 LMT  = LMT + 8
      IF (LMT.LT.LQMTE)            GO TO 22
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
