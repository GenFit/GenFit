*
* $Id: eropts.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: eropts.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.2  1996/07/25 15:01:02  ravndal
* CMZ to CVS conversion artefact removed
*
* Revision 1.1  1996/03/06 15:32:33  mclareni
* Add geane321 includes
*
*
#ifndef CERNLIB_ERPRCDES_EROPTS_INC
#define CERNLIB_ERPRCDES_EROPTS_INC
*
*      Sequence for conditions to stop tracking
*
* eropts.inc
*
      CHARACTER*8     CHOPTI
      LOGICAL         LEEXAC, LELENG, LEONLY, LEPLAN, LEPOIN, LEVOLU
      REAL            ERPLI, ERPLO, ERLENG
      INTEGER         NAMEER, NUMVER, IOVLER
      COMMON /EROPTS/ ERPLI(3,2), ERPLO(3,4,MXPRED), ERLENG(MXPRED),
     +                NAMEER(MXPRED), NUMVER(MXPRED), IOVLER(MXPRED),
     +                LEEXAC, LELENG, LEONLY, LEPLAN, LEPOIN, LEVOLU
      COMMON /EROPTC/CHOPTI
*
* eropts.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant231/incbeg.inc"
     +    EROPTS
#include "geant321/incend.inc"
#endif

#endif
