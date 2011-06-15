*
* $Id: gcshpt.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcshpt.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:35  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCSHPT_INC
#define CERNLIB_GEANT321_GCSHPT_INC
*
*
* gcshpt.inc
*
#if !defined(CERNLIB_INCLUDE)
C  Shells are numbered from 1 to 24.
C  Shells used:
C               K,L1,L2,L3,M1,M2,M3,M4,M5
C               N1,N2,N3,N4,N5,N6,N7,
C               O1,O2,O3,O4,O5,P1,P2,P3
C   VARIABLES:
C     NSHLST - value of Z for which the shells starts to be present
C     N1ST   - pointer to K shell of a given Z (in ESHELL array)
C     NSHLLS - Number of used shells for a given Z
C     ESHELL - Shells potentials in eV !!!
      INTEGER LENGTH,MAXSHL
      PARAMETER (LENGTH=  1409)
      PARAMETER (MAXSHL=24)
      INTEGER NSHLST,N1ST,NSHLLS
      REAL ESHELL
      DIMENSION NSHLST(MAXSHL),N1ST(MAXELZ),NSHLLS(MAXELZ)
      DIMENSION ESHELL(LENGTH)
      COMMON /GCSHPT/NSHLST,N1ST,NSHLLS,ESHELL
C
#endif
*
* gcshpt.inc
*
#if defined(CERNLIB_INCLUDE)
*
#include "geant321/incbeg.inc"
     +   gcshpt
#include "geant321/incend.inc"
#endif

#endif
