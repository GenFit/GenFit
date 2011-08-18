*
* $Id: gpdk.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gpdk.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:18  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      FUNCTION GPDK(A,B,C)
C
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C CALLED BY : PHASP
C ORIGIN : H.FESEFELDT (27-OCT-1983)
C
C     GPDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
C
      A2 = A*A
      B2 = B*B
      C2 = C*C
      IF(A2) 21,21,61
   61 CONTINUE
      ARG=A2+(B2-C2)**2/A2-2.0*(B2+C2)
      IF (ARG) 21,21,31
   21 GPDK=0.0
      GOTO 41
   31 CONTINUE
      GPDK = 0.5*SQRT(ABS(A2 + (B2-C2)**2/A2 - 2.0*(B2+C2)))
   41 CONTINUE
      RETURN
      END
