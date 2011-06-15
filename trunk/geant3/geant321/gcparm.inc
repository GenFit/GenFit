*
* $Id: gcparm.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcparm.inc,v $
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
#ifndef CERNLIB_GEANT321_GCPARM_INC
#define CERNLIB_GEANT321_GCPARM_INC
*
*
* gcparm.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCPARM/IPARAM,PCUTGA,PCUTEL,PCUTNE,PCUTHA,PCUTMU
     +             ,NSPARA,MPSTAK,NPGENE
      REAL PACUTS(5)
      EQUIVALENCE (PACUTS(1),PCUTGA)
#if defined(CERNLIB_TYPE)
#include "geant321/gtparm.inc"
#endif
C
#endif
*
* gcparm.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcparm
#include "geant321/incend.inc"
#endif

#endif
