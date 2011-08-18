*
* $Id: gcmate.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcmate.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:32  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCMATE_INC
#define CERNLIB_GEANT321_GCMATE_INC
*
*
* gcmate.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCMATE/NMAT,NAMATE(5),A,Z,DENS,RADL,ABSL
C
      INTEGER NMAT,NAMATE
      REAL A,Z,DENS,RADL,ABSL
C
#endif
*
* gcmate.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcmate
#include "geant321/incend.inc"
#endif

#endif
