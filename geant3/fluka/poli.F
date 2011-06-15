*
* $Id: poli.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: poli.F,v $
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
*$ CREATE POLI.FOR
*COPY POLI
*
*=== poli =============================================================*
*
      SUBROUTINE POLI(CS,SI)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      REAL RNDM(1)
C
C--------------------------------------------------
C*** RANDOM CHOICE OF ANGLE TETA (CS = COS(TETA),SI = SIN(TETA)
C--------------------------------------------------
      CALL GRNDM(RNDM,1)
      CS = 2.D+00 * RNDM(1) - 1.D+00
      SI = SQRT( ( 1.D+00 - CS ) * ( 1.D+00 + CS ) )
      RETURN
      END
