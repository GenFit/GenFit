*
* $Id: sverl2.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: sverl2.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:19  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.40  by  S.Giani
*-- Author :
      SUBROUTINE SVERL2(M,X)
      DIMENSION U(2)
      CALL GRNDM(U,2)
C PH      X=0.62666*LOG((1.+U(1))/(1.-U(1)))
      X=0.62666*LOG( (2.-U(1))/U(1) )
      IF(U(2).LT.0.5) X=-X
      X=M+X*SQRT(M*1.)
      RETURN
      END
