*
* $Id: gccuts.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gccuts.inc,v $
* Revision 1.2  2007/03/31 13:48:01  brun
*  From Andrea Fontana and Alberto Rotondi
* Added one extra parameter (GCALPHA) to the COMMON GCCUTS to drive the
* Urban/Landau models in the ERLAND routine.
*
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
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
#ifndef CERNLIB_GEANT321_GCCUTS_INC
#define CERNLIB_GEANT321_GCCUTS_INC
*
*
* gccuts.inc
*
* modified by A.Rotondi & A.Fontana 
* march 2007
* (reintroduced GCUTS(5))
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCCUTS/CUTGAM,CUTELE,CUTNEU,CUTHAD,CUTMUO,BCUTE,BCUTM
     +             ,DCUTE ,DCUTM ,PPCUTM,TOFMAX,GCUTS(5)
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtcuts.inc"

#endif
#endif
*
* gccuts.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gccuts
#include "geant321/incend.inc"
#endif

#endif
