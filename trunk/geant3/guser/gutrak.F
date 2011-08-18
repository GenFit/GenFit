*
* $Id: gutrak.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gutrak.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
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
      SUBROUTINE GUTRAK
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to control tracking of one track            *
C.    *                                                                *
C.    *    ==>Called by : GTREVE                                       *
C.    *                                                                *
C.    ******************************************************************
C.
C.
C.    ------------------------------------------------------------------
C.
*
      CALL G3TRACK
*                                                             END GUTRAK
      END
