*
* $Id: mziocf.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mziocf.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:26  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:47:21  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE MZIOCF (JUP,MXVAL)

C-    Find largest possible size of next group, slave to MZIOCH

#include "zebra/zkrakc.inc"
#include "zebra/quest.inc"
C--------------    END CDE                             --------------
      DIMENSION    MU(99)
      EQUIVALENCE (MU(1),IQHOLK(1))

      EQUIVALENCE (NGR,IQUEST(1)), (NGRU,IQUEST(2))

      DIMENSION    JUP(9), MXVAL(9)


      JU  = JUP(1)
      MXC = MU(JU+2)

      DO 24  JL=2,NGR
      JU = JU + 2
      MXC = MAX (MU(JU+2),MXC)
      IF (MXC.GE.MXVAL(JL))        GO TO 29
   24 CONTINUE
      NGRU = NGR
      RETURN

   29 NGRU = JL - 1
      RETURN
      END
*      ==================================================
#include "zebra/qcardl.inc"
