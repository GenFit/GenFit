*
* $Id: gcflag.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcflag.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:31  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCFLAG_INC
#define CERNLIB_GEANT321_GCFLAG_INC
*
*
* gcflag.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtflag.inc"
#endif
#endif
*
* gcflag.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcflag
#include "geant321/incend.inc"
#endif

#endif
