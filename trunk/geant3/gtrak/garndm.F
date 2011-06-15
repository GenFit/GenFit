*
* $Id: garndm.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: garndm.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:40  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      FUNCTION GARNDM(DUMMY)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       GARNDM = -LOG(RNDM(DUMMY))                               *
C.    *                                                                *
C.    *    ==>Called by : GBREME,GBREMM,GCOMP,GDRAY,GHADRO,GHEISH,     *
C.    *                   GLANDO,GPAIRG,GPAIRM,GRAYL,GLTRAC            *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION RNDM(1)
C.
C.    ------------------------------------------------------------------
C.
      CALL GRNDM(RNDM,1)
      GARNDM=-LOG(RNDM(1))
      END
