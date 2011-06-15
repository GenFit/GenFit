*
* $Id: altra.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: altra.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:36  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE ALTRA.FOR
*COPY ALTRA
*
*=== altra ============================================================*
*
      SUBROUTINE ALTRA(GA,BGX,BGY,BGZ,PCX,PCY,PCZ,EC,P,PX,PY,PZ,E)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C*** ARBITRARY LORENTZ TRANSFORM
      EP = PCX * BGX + PCY * BGY + PCZ * BGZ
      E  = GA * EC + EP
      PE = EP / (GA + 1.D0) + EC
      PX = PCX + BGX * PE
      PY = PCY + BGY * PE
      PZ = PCZ + BGZ * PE
      TMP30 = 1.D30
      PX = MIN(TMP30, MAX(-TMP30,PX))
      PY = MIN(TMP30, MAX(-TMP30,PY))
      PZ = MIN(TMP30, MAX(-TMP30,PZ))
      P  = SQRT (PX * PX + PY * PY + PZ * PZ)
      RETURN
      END
