*
* $Id: guskip.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: guskip.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUSKIP (ISKIP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   User routine to skip unwanted tracks                         *
C.    *                                                                *
C.    *   Called by : GSSTAK                                           *
C.    *   Author    : F.Bruyant                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
C.    ------------------------------------------------------------------
*
      ISKIP = 0
*                                                             END GUSKIP
      END
 
