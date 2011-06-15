*
* $Id: sbyt.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: sbyt.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:14  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/sbyt.F"
#else
      SUBROUTINE SBYT (IT,IZW,IZP,NZB)
C
C CERN PROGLIB# M421    SBYT            .VERSION KERNFOR  4.23  891215
C MOD. true default 24/2/89, JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine
C
#include "kerngen/wordsize.inc"
#include "kerngen/iall1s.inc"
#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

#if !defined(CERNLIB_QISASTD)
      MSK = ISHFTR (IALL11, NBITPW-NZB)
      IZW = IOR ( IAND (IZW, NOT(ISHFTL(MSK,IZP-1)))
     +,           ISHFTL(IAND(IT,MSK),IZP-1))
#endif
#if defined(CERNLIB_QISASTD)
      MSK = ISHFT  (IALL11,-NBITPW+NZB)
      IZW = IOR ( IAND (IZW, NOT(ISHFT (MSK,IZP-1)))
     +,           ISHFT (IAND(IT,MSK),IZP-1))
#endif
      RETURN
      END
#endif
