*
* $Id: erwork.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: erwork.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.2  1996/07/25 15:01:16  ravndal
* CMZ to CVS conversion artefact removed
*
* Revision 1.1  1996/03/06 15:32:34  mclareni
* Add geane321 includes
*
*
#ifndef CERNLIB_ERPRCDES_ERWORK_INC
#define CERNLIB_ERPRCDES_ERWORK_INC
*
*      Working sequences for the interface
*
* erwork.inc
*
      DOUBLE PRECISION EI, EF, ASDSC
      COMMON /ERWORK/ EI(15), EF(15), ASDSC(5,5),
     +                   XI(3), PPI(3), HI(9),
     +                   XF(3), PF(3),  HF(9),
     +                   CHTR, DEDX2, BACKTR, CUTEK, TLGCM2, TLRAD
*
* erwork.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
     +    ERWORK
#include "geant321/incend.inc"
#endif

#endif
