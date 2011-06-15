*
* $Id: guplsh.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: guplsh.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:47  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      FUNCTION GUPLSH(MED0,MED1)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *    ==>Called by : GLISUR                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctmed.inc"
C.
C.    ------------------------------------------------------------------
C.
*
* *** By default this defines perfect smoothness
      GUPLSH = 1.
C
      END
