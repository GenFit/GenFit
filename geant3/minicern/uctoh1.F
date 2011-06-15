*
* $Id: uctoh1.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: uctoh1.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.2  1996/04/02 23:17:46  thakulin
* Add support for EPC Fortran:  remove char-int equivalences and use
* F90 transfer facility instead.
*
* Revision 1.1.1.1  1996/02/15 17:50:15  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UCTOH1 (MS,MT,NCH)
C
C CERN PROGLIB# M409    UCTOH1          .VERSION KERNFOR  4.21  890323
C ORIG. 10/02/89  JZ
C
      CHARACTER    MS*99
      DIMENSION    MT(99)

#include "kerngen/wordsize.inc"
      CHARACTER    CHWORD*(NCHAPW)
      CHARACTER    BLANK *(NCHAPW)
      PARAMETER   (BLANK = ' ')
      INTEGER      IWORD
#if !defined(CERNLIB_F90) && !defined(CERNLIB_QFEPC)
      EQUIVALENCE (IWORD,CHWORD)
#endif

      IF   (NCH)             91, 29, 11

   11 CHWORD = BLANK

      DO 24 JS=1,NCH
      CHWORD(1:1) = MS(JS:JS)
#if defined(CERNLIB_F90) || defined(CERNLIB_QFEPC)
      MT(JS)      = transfer(CHWORD,IWORD)
#else
      MT(JS)      = IWORD
#endif
   24 CONTINUE
   29 RETURN

   91 CALL ABEND
      END
