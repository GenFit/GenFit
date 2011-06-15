*
* $Id: frhinc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: frhinc.F,v $
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
*$ CREATE FRHINC.FOR
*COPY FRHINC
*
*=== frhinc ===========================================================*
*
      FUNCTION FRHINC (R)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
*
#include "geant321/nucstf.inc"
*
      IF ( R .LE. RADIU0 ) THEN
         FRHINC = RICORE (R)
      ELSE IF ( R .LE. RADIU1 ) THEN
         FRHINC = RISKIN (R)
      ELSE IF ( R .LE. RADTOT ) THEN
         FRHINC = RIHALO (R)
      ELSE
         FRHINC = RIHALO (RADTOT)
      END IF
      RETURN
      END
