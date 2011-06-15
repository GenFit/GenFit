*
* $Id: fradnc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: fradnc.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:43  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:22  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:22:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE FRADNC.FOR
*COPY FRADNC
*
*=== fradnc ===========================================================*
*
      FUNCTION FRADNC (RHO)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
*
#include "geant321/nucstf.inc"
*
      IF ( RHO .GE. RHOCOR ) THEN
         FRADNC = RACORE (RHO)
      ELSE IF ( RHO .GE. RHOSKN ) THEN
         FRADNC = RASKIN (RHO)
      ELSE IF ( RHO .GT. 0.D+00 ) THEN
         FRADNC = RAHALO (RHO)
      ELSE
         FRADNC = RADTOT
      END IF
      RETURN
      END
