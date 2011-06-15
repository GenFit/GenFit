*
* $Id: trcom3.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: trcom3.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.2  1996/07/25 15:01:26  ravndal
* CMZ to CVS conversion artefact removed
*
* Revision 1.1  1996/03/06 15:32:35  mclareni
* Add geane321 includes
*
*
#ifndef CERNLIB_ERPRCDES_TRCOM3_INC
#define CERNLIB_ERPRCDES_TRCOM3_INC
*
*      Sequences of the EMC routines
*
* trcom3.inc
*
      DOUBLE PRECISION A, B, S, TN, T, COSL, SINL, COSP, SINP, COSL1
      INTEGER          NEW
      COMMON /TRCOM3/ A(5,5), B(5,5), S(15),
     +                TN(3), T(5), COSL, SINL, COSP, SINP, COSL1, NEW
*
* trcom3.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
     +    TRCOM3
#include "geant321/incend.inc"
#endif

#endif
