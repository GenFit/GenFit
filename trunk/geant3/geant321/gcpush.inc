*
* $Id: gcpush.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcpush.inc,v $
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
#ifndef CERNLIB_GEANT321_GCPUSH_INC
#define CERNLIB_GEANT321_GCPUSH_INC
*
*
* gcpush.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCPUSH/NCVERT,NCKINE,NCJXYZ,NPVERT,NPKINE,NPJXYZ
      INTEGER       NCVERT,NCKINE,NCJXYZ,NPVERT,NPKINE,NPJXYZ
C
#endif
*
* gcpush.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcpush
#include "geant321/incend.inc"
#endif

#endif
