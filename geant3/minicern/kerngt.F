*
* $Id: kerngt.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: kerngt.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.2  1997/02/04 17:36:16  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:31:21  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:40  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_WINNT)
#include "wntgs/kerngt.F"
#elif defined(CERNLIB_QMHPX)
#include "hpxgs/kerngt.F"
#elif defined(CERNLIB_QMIRT)||defined(CERNLIB_QMIRTD)
#include "irtgs/kerngt.F"
#elif defined(CERNLIB_QMLNX)
#include "lnxgs/kerngt.F"
#elif defined(CERNLIB_QMSGI)
#include "sgigs/kerngt.F"
#elif defined(CERNLIB_QMSUN)
#include "sungs/kerngt.F"
#elif defined(CERNLIB_QMVMI)||defined(CERNLIB_QMVAO)
#include "allgs/kerngt.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/kerngt.F"
#elif defined(CERNLIB_QMDOS)
#include "dosgs/kerngt.F"
#elif defined(CERNLIB_QMMPW)
#include "mpwgs/kerngt.F"
#elif defined(CERNLIB_QMOS9)
#include "os9gs/kerngt.F"
#else
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNFOR  4.29  910718
C ORIG. 18/12/78
C
      COMMON /SLATE/ ICH(32),DUM(8)

      CALL VBLANK (ICH,40)
      CALL GSHOWF (LUN)
      RETURN
      END
#endif
