*
* $Id: gcvolu.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcvolu.inc,v $
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
#ifndef CERNLIB_GEANT321_GCVOLU_INC
#define CERNLIB_GEANT321_GCVOLU_INC
*
*
* gcvolu.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15),
     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15),
     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3)
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtvolu.inc"
#endif
#endif
*
* gcvolu.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcvolu
#include "geant321/incend.inc"
#endif

#endif
