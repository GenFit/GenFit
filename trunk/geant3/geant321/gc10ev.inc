*
* $Id: gc10ev.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gc10ev.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GC10EV_INC
#define CERNLIB_GEANT321_GC10EV_INC
*
*
* gc10ev.inc
*
#if !defined(CERNLIB_INCLUDE)
      REAL G10EV,TENEV
      PARAMETER (G10EV=1.0E-8)
      PARAMETER (TENEV=1.E-2)
C
#endif
*
* gc10ev.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gc10ev
#include "geant321/incend.inc"
#endif

#endif
