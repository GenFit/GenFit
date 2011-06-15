*
* $Id: cmass.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cmass.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:44  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CMASS_INC
#define CERNLIB_GEANT321_CMASS_INC
*
*
* cmass.inc
*
#include "geant321/camass.inc"
      DIMENSION PMASS(12)
      EQUIVALENCE(PMASS(1),XMASS(0))

#endif
