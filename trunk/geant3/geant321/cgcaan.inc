*
* $Id: cgcaan.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cgcaan.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGCAAN_INC
#define CERNLIB_GEANT321_CGCAAN_INC
*
*
* cgcaan.inc
*
#if !defined(CERNLIB_INCLUDE)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION RAD,COSI,SINI,COSDA(2),SINDA(2),XYHA(4,2)
      INTEGER          NASTP(2)
      COMMON /CGCAAN/  RAD,COSI,SINI,COSDA,SINDA,XYHA,NASTP
     +                ,IFULL,ATREDG
#endif
#endif
#if !defined(CERNLIB_INCLUDE)
#if defined(CERNLIB_SINGLE)
      REAL             COSDA(2),SINDA(2),XYHA(4,2)
      INTEGER          NASTP(2)
      COMMON /CGCAAN/  RAD,COSI,SINI,COSDA,SINDA,XYHA,NASTP
     +                ,IFULL,ATREDG
#endif
#endif
*
* cgcaan.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGCAAN
#include "geant321/incend.inc"
*
#endif

#endif
