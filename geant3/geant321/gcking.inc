*
* $Id: gcking.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcking.inc,v $
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
#ifndef CERNLIB_GEANT321_GCKING_INC
#define CERNLIB_GEANT321_GCKING_INC
*
*
* gcking.inc
*
#if !defined(CERNLIB_INCLUDE)
#include "geant321/gckmax.inc"
      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN),
     +                           TOFD(MXGKIN),IFLGK(MXGKIN)
      INTEGER       KCASE,NGKINE ,IFLGK,MXPHOT,NGPHOT
      REAL          GKIN,TOFD,XPHOT
C
      PARAMETER (MXPHOT=800)
      COMMON/GCKIN2/NGPHOT,XPHOT(11,MXPHOT)
C
      COMMON/GCKIN3/GPOS(3,MXGKIN)
      REAL          GPOS
C
#endif
*
* gcking.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcking
#include "geant321/incend.inc"
#endif

#endif
