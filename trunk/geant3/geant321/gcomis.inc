*
* $Id: gcomis.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcomis.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:34  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCOMIS_INC
#define CERNLIB_GEANT321_GCOMIS_INC
*
*
* gcomis.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCOMIS/JUINIT,JUGEOM,JUKINE,JUSTEP,JUOUT,JULAST
      DIMENSION JPCOMS(6)
      EQUIVALENCE (JPCOMS,JUINIT)
*
#endif
*
* gcomis.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcomis
#include "geant321/incend.inc"
#endif

#endif
