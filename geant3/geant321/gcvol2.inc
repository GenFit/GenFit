*
* $Id: gcvol2.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcvol2.inc,v $
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
#ifndef CERNLIB_GEANT321_GCVOL2_INC
#define CERNLIB_GEANT321_GCVOL2_INC
*
*
* gcvol2.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCVOL2/NLEVE2,NAMES2(15),NUMB2(15),
     +LVOL2(15),LIND2(15),INFRO2,NLDEV2(15),LINMX2(15),
     +GTRAN2(3,15),GRMAT2(10,15),GONLY2(15),GLX2(15)
      INTEGER NLEVE2,NAMES2,NUMB2,LVOL2,LIND2,INFRO2,NLDEV2,LINMX2
      REAL GTRAN2,GRMAT2,GONLY2,GLX2
C
#endif
*
* gcvol2.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcvol2
#include "geant321/incend.inc"
#endif

#endif
