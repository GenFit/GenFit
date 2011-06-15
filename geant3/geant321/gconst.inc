*
* $Id: gconst.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gconst.inc,v $
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
#ifndef CERNLIB_GEANT321_GCONST_INC
#define CERNLIB_GEANT321_GCONST_INC
*
*
* gconst.inc
*
#if !defined(CERNLIB_INCLUDE)
#if defined(CERNLIB_TYPE)
#include "geant321/gtcons.inc"
#endif
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
#endif
*
* gconst.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gconst
#include "geant321/incend.inc"
#endif

#endif
