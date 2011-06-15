*
* $Id: gcjump.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcjump.inc,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
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
#ifndef CERNLIB_GEANT321_GCJUMP_INC
#define CERNLIB_GEANT321_GCJUMP_INC
*
*
* gcjump.inc
*
#if !defined(CERNLIB_INCLUDE)
#if defined(CERNLIB_TYPE)
#include "geant321/gtjump.inc"
#endif
      PARAMETER    (MAXJMP=30)
#if defined(CERNLIB_LXIA64) && !defined(CERNLIB_TYPE)
      INTEGER*8     JUDCAY, JUDIGI, JUDTIM, JUFLD , JUHADR, JUIGET,
     +              JUINME, JUINTI, JUKINE, JUNEAR, JUOUT , JUPHAD,
     +              JUSKIP, JUSTEP, JUSWIM, JUTRAK, JUTREV, JUVIEW,
     +              JUPARA
      INTEGER*8 JMPADR
#endif

      COMMON/GCJUMP/JUDCAY, JUDIGI, JUDTIM, JUFLD , JUHADR, JUIGET,
     +              JUINME, JUINTI, JUKINE, JUNEAR, JUOUT , JUPHAD,
     +              JUSKIP, JUSTEP, JUSWIM, JUTRAK, JUTREV, JUVIEW,
     +              JUPARA
      DIMENSION     JMPADR(MAXJMP)
      EQUIVALENCE  (JMPADR(1), JUDCAY)
*
#endif
*
* gcjump.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcjump
#include "geant321/incend.inc"
#endif

#endif
