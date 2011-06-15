*
* $Id: gctime.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gctime.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCTIME_INC
#define CERNLIB_GEANT321_GCTIME_INC
*
*
* gctime.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCTIME/TIMINT,TIMEND,ITIME,IGDATE,IGTIME
      INTEGER ITIME,IGDATE,IGTIME
      REAL TIMINT,TIMEND
C
#endif
*
* gctime.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gctime
#include "geant321/incend.inc"
#endif

#endif
