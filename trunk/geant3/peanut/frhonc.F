*
* $Id: frhonc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: frhonc.F,v $
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
*$ CREATE FRHONC.FOR
*COPY FRHONC
*
*=== frhonc ===========================================================*
*
      FUNCTION FRHONC (R)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
*
#include "geant321/nucstf.inc"
*
      IF ( R .LE. RADIU0 ) THEN
         FRHONC = RHCORE (R)
      ELSE IF ( R .LE. RADIU1 ) THEN
         FRHONC = RHSKIN (R)
      ELSE IF ( R .LE. RADTOT ) THEN
         FRHONC = RHHALO (R)
      ELSE
         FRHONC = 0.D+00
      END IF
      RETURN
      END
