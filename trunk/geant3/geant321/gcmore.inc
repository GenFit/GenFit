*
#ifndef CERNLIB_GEANT321_GCMORE_INC
#define CERNLIB_GEANT321_GCMORE_INC
*
*
* gcmore.inc
*
* introduced by A.Rotondi & A.Fontana 
* may 2007
* 
#if !defined(CERNLIB_INCLUDE)
      COMMON/GCMORE/GCALPHA,ICLOSE,PFINAL(3),DSTRT,WIRE1(3),WIRE2(3),
     +              P1(3),P2(3),P3(3),CLENG(3)
C
C     input to ERLAND:
C     GCALPHA: energy cut parameter for energy loss fluctuations
C
C     input to EUSTEP:
C     ICLOSE: =1 the use of the common is enabled for a point
C             =2 the use of the common is enabled for a wire
C             =0 the common is empty
C     PFINAL(3): assigned point
C     WIRE1(3): first point of a wire
C     WIRE2(3): second point of a wire
C     DSTRT: assigned distance between initial point in ERTRAK 
C            and PFINAL
C
C     output written by EUSTEP:
C     P1(3): point previous to the point of closest approach to 
C            PFINAL() or wire
C     P2(3): point of closest approach to PFINAL() or wire
C     P3(3): point next to the point of closest approach to 
C            PFINAL() or wire
C     CLENG(3): track length to the previous 3 points
C 
#if defined(CERNLIB_TYPE)
#include "geant321/gtmore.inc"

#endif
#endif
*
* gcmore.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcmore
#include "geant321/incend.inc"
#endif

#endif
