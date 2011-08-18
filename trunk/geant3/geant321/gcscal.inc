*
* $Id: gcscal.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcscal.inc,v $
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
#ifndef CERNLIB_GEANT321_GCSCAL_INC
#define CERNLIB_GEANT321_GCSCAL_INC
*
*
* gcscal.inc
*
#if !defined(CERNLIB_INCLUDE)
#if defined(CERNLIB_TYPE)
#include "geant321/gtscal.inc"
#endif
      PARAMETER(MXSLNK=100)
      COMMON/GCSCAL/ ISLINK(MXSLNK)
      EQUIVALENCE (LSLAST,ISLINK(MXSLNK))
      EQUIVALENCE (LSCAN ,ISLINK(1)),(LSTEMP,ISLINK(2))
      EQUIVALENCE (LSPARA,ISLINK(3)),(LSERAY,ISLINK(4))
*
#endif
*
* gcscal.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcscal
#include "geant321/incend.inc"
#endif

#endif
