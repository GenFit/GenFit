*
* $Id: cterp.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cterp.F,v $
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
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE CTERP(X1,X2,X,Y1,Y2,Y)
C       THIS ROUTINE PERFORMS LINEAR INTERPOLATION
      Y=Y2-(X2-X)*(Y2-Y1)/(X2-X1)
      RETURN
      END
