*
* $Id: svgam7.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: svgam7.F,v $
* Revision 1.2  2005/08/11 13:34:26  brun
* From Peter Hristov:
* an additional protection to avoid random numbers exactly equal to 1 or 0.
*
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:19  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:12  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 19/05/94  14.28.55  by  S.Ravndal
*-- Author :
      SUBROUTINE SVGAM7(A,X)
      DIMENSION U(2)
      REAL LA
      GA=A-1.
      LA=SQRT(2.*A-1.)
      EP=1.570796327+ATAN(GA/LA)
      RO=1.570796327-EP
    1 CALL GRNDM(U,2)
      IF ((U(1) .GE. 1.0) .OR. (U(2) .LE. 0.0) ) GOTO 1
      X=GA+LA*TAN(EP*U(1)+RO)
      IF (X.LE.0.) GOTO 1
      Y=LOG(1.+((X-GA)/LA)**2) +GA*LOG(X/GA)-X+GA
      IF(LOG(U(2)).GT.Y) GOTO 1
      RETURN
      END
