*
* $Id: trafo.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: trafo.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE TRAFO.FOR
*COPY TRAFO
*
*=== trafo ============================================================*
*
      SUBROUTINE TRAFO(GAM,BGAM,CX,CY,CZ,COD,COF,SIF,P,ECM,
     1PL,CXL,CYL,CZL,EL)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C     LORENTZ TRANSFORMATION INTO THE LAB - SYSTEM
      SID=SQRT(1.D0-COD*COD)
      PLX=P*SID*COF
      PLY=P*SID*SIF
      PCMZ=P*COD
      PLZ=GAM*PCMZ+BGAM*ECM
      PL=SQRT(PLX*PLX+PLY*PLY+PLZ*PLZ)
      EL=GAM*ECM+BGAM*PCMZ
C     ROTATION INTO THE ORIGINAL DIRECTION
      COZ=PLZ/PL
      SIZ=SQRT(1.D0-COZ**2)
      CALL TTRANS(CX,CY,CZ,COZ,SIZ,SIF,COF,CXL,CYL,CZL)
      RETURN
      END
