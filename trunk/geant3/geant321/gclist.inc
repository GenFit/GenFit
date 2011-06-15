*
* $Id: gclist.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gclist.inc,v $
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
#ifndef CERNLIB_GEANT321_GCLIST_INC
#define CERNLIB_GEANT321_GCLIST_INC
*
*
* gclist.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
#if defined(CERNLIB_TYPE)
#include "geant321/gtlist.inc"
#endif
#endif
*
* gclist.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gclist
#include "geant321/incend.inc"
#endif

#endif
