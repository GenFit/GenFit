*
* $Id: nfac.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: nfac.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:18  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      INTEGER FUNCTION NFAC(N)
C
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C ORIGIN : H.FESEFELDT (27-OCT-1983)
C
      NFAC=1.
      M=N
      IF(M.LE.1) RETURN
      IF(M.GT.10) M=10
      DO 1 I=2,M
    1 NFAC=NFAC*I
      RETURN
      END
