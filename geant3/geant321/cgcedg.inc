*
* $Id: cgcedg.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cgcedg.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGCEDG_INC
#define CERNLIB_GEANT321_CGCEDG_INC
*
*
* cgcedg.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (INPART = 30, NTTT = 50)
      COMMON /CGCEDG/ TDEL,NPART,TSTRT(INPART),TEND(INPART),
     +                KX,KY,KZ,NPART2,TSTRT2(INPART),TEND2(INPART),
     +                XA,YA,ZA,XB,YB,ZB,XDELT,YDELT,ZDELT,
     +                ABC(3,2),DEDGE(5),TTT(2,NTTT),ITTT(NTTT)
      REAL AA(3),BB(3),DELTA(3)
      EQUIVALENCE     (AA,XA),(BB,XB),(DELTA,XDELT)
#endif
*
* cgcedg.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGCEDG.
#include "geant321/incend.inc"
*
#endif

#endif
