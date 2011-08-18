*
* $Id: beta.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: beta.F,v $
* Revision 1.2  2003/06/10 16:12:15  brun
* From Andreas Morsch:
* I changed the name of the function BETA from geant3/fluka/beta.F to
* BETA_G3. The reason is that if one uses PDF and Geant3 in the same
* simulation there is a clash between a similar function in PDF.
*
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:36  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE BETA.FOR
*COPY BETA
*
*=== beta =============================================================*
*
      FUNCTION BETA_G3(X1,X2,BET)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      AX=0.D0
      BETX1=BET*X1
      IF(BETX1.LT.70.D0) AX=-1.D0/BET**2*(BETX1+1.D0)*EXP(-BETX1)
      AY=1.D0/BET**2*(BET*X2+1.D0)*EXP(-BET*X2)
      BETA_G3=AX+AY
      RETURN
      END
