*
* $Id: glight.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: glight.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:04  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
*
      SUBROUTINE GLIGHT(XXLL,YYLL,ZZLL,IILL)
*
#include "geant321/gcrayt.inc"
*
      IOMBRA=10
      XLPOS=XXLL
      YLPOS=YYLL
      ZLPOS=ZZLL
      INTEN=IILL
      IF(INTEN.EQ.0)IOMBRA=0
*
      END
