*
* $Id: ranres.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ranres.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:19  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      FUNCTION RANRES(X)
C
C *** RESTRICTED RANDOM NUMBERS TO BE USED AS ARGUMENT IN LOG ETC... ***
C *** NVE 13-JUL-1988 CERN GENEVA ***
C
C NOTE : 0 < RANRES < 1
      DIMENSION RNDM(1)
C
 1    CONTINUE
      CALL GRNDM(RNDM,1)
      RANRES=RNDM(1)
      IF ((RANRES .LE. 0.) .OR. (RANRES .GE. 1.)) GO TO 1
      RETURN
      END
