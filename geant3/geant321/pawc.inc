*
* $Id: pawc.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: pawc.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_PAWC_INC
#define CERNLIB_GEANT321_PAWC_INC
*
*
* pawc.inc
*
#if !defined(CERNLIB_INCLUDE)
#if defined(CERNLIB_TYPE)
#include "geant321/pawct.inc"
#endif
      COMMON /PAWC/ NWPAW,IXPAWC,IHBOOK,IXHIGZ,IXKUIP,IFENCE(5),
     +              LLMAIN, WWS(9989)
      DIMENSION IQQ(2),QQ(2),LQQ(8000)
      EQUIVALENCE (QQ(1),IQQ(1),LQQ(9)),(LQQ(1),LLMAIN)
C
#endif
*
* pawc.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   pawc
#include "geant321/incend.inc"
#endif

#endif
