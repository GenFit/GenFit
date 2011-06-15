*
* $Id: gctmed.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gctmed.inc,v $
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
#ifndef CERNLIB_GEANT321_GCTMED_INC
#define CERNLIB_GEANT321_GCTMED_INC
*
*
* gctmed.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX
     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD
      COMMON/GCTLIT/THRIND,PMIN,DP,DNDL,JMIN,ITCKOV,IMCKOV,NPCKOV
C
#if defined(CERNLIB_TYPE)
#include "geant321/gttmed.inc"
#endif
#endif
*
* gctmed.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gctmed
#include "geant321/incend.inc"
#endif

#endif
