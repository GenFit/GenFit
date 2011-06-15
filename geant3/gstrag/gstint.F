*
* $Id: gstint.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gstint.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      FUNCTION G3STINT(X,A,F,N)
*
*       Simple linear interpolation
*
      DIMENSION A(N),F(N)
      DATA K1,K2,NO/0,0,0/
*
*       Check (N==NO) & (A(K1)<X<A(K2))
*
      IF (N.NE.NO.OR.A(K1).GT.X.OR.X.GT.A(K2)) THEN
*
*       Skip Unless N>1
*
         IF (N.LE.1) THEN
            G3STINT=F(1)
            GO TO 999
         ENDIF
         K1=1
         K2=N
         NO=N
*
*       Loop over K2-K1>1
*
  10     K=(K1+K2)/2
*
*       Check A(K)<X
*
         IF (A(K).LT.X) THEN
            K1=K
         ELSE
            K2=K
         END IF
         IF(K2-K1.GT.1) GO TO 10
      END IF
      X1=A(K1)
      X2=A(K2)
      G3STINT=(F(K1)*(X-X2)+F(K2)*(X1-X))/(X1-X2)
*
 999  END
