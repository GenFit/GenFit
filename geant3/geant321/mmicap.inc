*
* $Id: mmicap.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mmicap.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:45  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_MMICAP_INC
#define CERNLIB_GEANT321_MMICAP_INC
*
*
* mmicap.inc
*
C
#include "geant321/mmicab.inc"
#include "geant321/gcbank.inc"
      DIMENSION LD(1),D(1)
      EQUIVALENCE (D(1),Q(1)),(LD(1),IQ(1))
      CHARACTER*24 DATSTR
      CHARACTER*80 COMMEN
      COMMON/ MICMAT / DATSTR,COMMEN,MATIDS(100,20,2)
C

#endif
