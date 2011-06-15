*
* $Id: gcgobj.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcgobj.inc,v $
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
#ifndef CERNLIB_GEANT321_GCGOBJ_INC
#define CERNLIB_GEANT321_GCGOBJ_INC
*
*
* gcgobj.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (NTRCG=1)
      PARAMETER (NWB=207,NWREV=100,NWS=1500)
      PARAMETER (C2TOC1=7.7, C3TOC1=2.,TVLIM=1296.)
      COMMON /GCGOBJ/IST,IFCG,ILCG,NTCUR,NFILT,NTNEX,KCGST
     +             ,NCGVOL,IVFUN,IVCLOS,IFACST,NCLAS1,NCLAS2,NCLAS3
      COMMON /CGBLIM/IHOLE,CGXMIN,CGXMAX,CGYMIN,CGYMAX,CGZMIN,CGZMAX
C
#endif
*
* gcgobj.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcgobj
#include "geant321/incend.inc"
#endif

#endif
