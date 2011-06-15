*
* $Id: gsecti.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gsecti.inc,v $
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
#ifndef CERNLIB_GEANT321_GSECTI_INC
#define CERNLIB_GEANT321_GSECTI_INC
*
*
* gsecti.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GSECTI/ AIEL(20),AIIN(20),AIFI(20),AICA(20),ALAM,K0FLAG
      INTEGER K0FLAG
      REAL AIEL,AIIN,AIFI,AICA,ALAM
C
#endif
*
* gsecti.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gsecti
#include "geant321/incend.inc"
#endif

#endif
