*
* $Id: gcrz.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcrz.inc,v $
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
#ifndef CERNLIB_GEANT321_GCRZ_INC
#define CERNLIB_GEANT321_GCRZ_INC
*
*
* gcrz.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCRZ1/NRECRZ,NRGET,NRSAVE,LRGET(20),LRSAVE(20)
      INTEGER      NRECRZ,NRGET,NRSAVE,LRGET    ,LRSAVE
      COMMON/GCRZ2/RZTAGS
      CHARACTER*8 RZTAGS(4)
C
#endif
*
* gcrz.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcrz
#include "geant321/incend.inc"
#endif

#endif
