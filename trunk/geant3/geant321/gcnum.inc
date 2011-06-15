*
* $Id: gcnum.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcnum.inc,v $
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
#ifndef CERNLIB_GEANT321_GCNUM_INC
#define CERNLIB_GEANT321_GCNUM_INC
*
*
* gcnum.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
     +            ,NSTMAX,NVERTX,NHEAD,NBIT
      COMMON /GCNUMX/ NALIVE,NTMSTO
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtnum.inc"
#endif
#endif
*
* gcnum.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcnum
#include "geant321/incend.inc"
#endif

#endif
