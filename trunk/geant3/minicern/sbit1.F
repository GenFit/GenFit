*
* $Id: sbit1.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: sbit1.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:19  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/sbit1.F"
#else
      SUBROUTINE SBIT1 (IZW,IZP)
C
C CERN PROGLIB# M421    SBIT1           .VERSION KERNFOR  4.23  891215
C MOD. true default 24/2/89, JZ
C
C     This non-ANSI code is a default which may be slow, if so
C     it should be replaced by a machine-specific fast routine

#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

#if !defined(CERNLIB_QISASTD)
      IZW = IOR (IZW, ISHFTL(1,IZP-1))
#endif
#if (defined(CERNLIB_QISASTD))&&(!defined(CERNLIB_QMILSTD))
      IZW = IOR (IZW, ISHFT (1,IZP-1))
#endif
#if defined(CERNLIB_QMILSTD)
      IZW = IBSET (IZW, IZP-1)
#endif
      RETURN
      END
#endif
