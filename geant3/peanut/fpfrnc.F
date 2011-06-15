*
* $Id: fpfrnc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: fpfrnc.F,v $
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
*$ CREATE FPFRNC.FOR
*COPY FPFRNC
*
*=== fpfrnc ===========================================================*
*
      FUNCTION FPFRNC (RHO,I)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
#include "geant321/nucgeo.inc"
*
#include "geant321/nucstf.inc"
*
      IF ( RHO .GE. RHOCOR ) THEN
         FPFRNC = PFCORE (RHO,I)
      ELSE IF ( RHO .GE. RHOSKN ) THEN
         FPFRNC = PFSKIN (RHO,I)
      ELSE IF ( RHO .GT. 0.D+00 ) THEN
         FPFRNC = PFHALO (RHO,I)
      ELSE
         FPFRNC = 0.D+00
      END IF
      RETURN
      END
