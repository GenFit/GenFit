*
* $Id: s_nucio.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: s_nucio.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:44  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_S_NUCIO_INC
#define CERNLIB_GEANT321_S_NUCIO_INC
*
*
* s_nucio.inc
*
      COMMON/NUCIN /TECM,AMASS(18),NPG,KGENEV
      COMMON/NUCOUT/PCM(5,18),WGT
#if !defined(CERNLIB_SINGLE)
#include "geant321/nuciod.inc"
#endif
C

#endif
