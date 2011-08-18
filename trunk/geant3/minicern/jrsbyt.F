*
* $Id: jrsbyt.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: jrsbyt.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:12  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/jrsbyt.F"
#else
      FUNCTION JRSBYT (VAL,X,JTH,NBITS)
C
C CERN PROGLIB# M421    JRSBYT          .VERSION KERNFOR  4.07  830624
C ORIG. 03/05/72 JZ
C
      INTEGER   VAL, X
C
#include "kerngen/q_jbit.inc"
C
      JRSBYT= JBYT (X,JTH,NBITS)
      CALL SBYT (VAL,X,JTH,NBITS)
      RETURN
      END
#endif
