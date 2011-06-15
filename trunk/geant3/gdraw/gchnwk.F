*
* $Id: gchnwk.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gchnwk.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:03  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
*
      SUBROUTINE GCHNWK
*
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcspee.inc"
*
      SAVE FIRST
      LOGICAL FIRST
      DATA FIRST/.TRUE./
*
      IF(FIRST)THEN
         FIRST=.FALSE.
         CALL IOPWK(2,1,9)
      ENDIF
      CALL IDAWK(1)
      CALL IACWK(2)
*
 999  END
