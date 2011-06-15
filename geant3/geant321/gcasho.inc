*
* $Id: gcasho.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcasho.inc,v $
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
#ifndef CERNLIB_GEANT321_GCASHO_INC
#define CERNLIB_GEANT321_GCASHO_INC
*
*
* gcasho.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCASHO/ZMED,AMED,DMED,E0MED,ZSMED(50),ESMED(50),ALFA,
     *             STEP,PLIN,PLOG,BE2,PLASM,TRNSMA,
     *             BOSC(50),AOSC(50),EOSC(50),ZOSC(50),EMEAN,
     *             CMGO(2000),EMGO,EMGOMI,
     *             NSMED,IOSC(50),NOSC,NMGO,NMGOMA
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtasho.inc"
#endif
#endif

#endif
