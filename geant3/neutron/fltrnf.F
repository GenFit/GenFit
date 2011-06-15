*
* $Id: fltrnf.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: fltrnf.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:43  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
C*********************************************************************
      FUNCTION FLTRNF(X)
C*********************************************************************
#include "geant321/crandm.inc"
      FLTRNF = RANDC(ISEED)
      RETURN
      END
