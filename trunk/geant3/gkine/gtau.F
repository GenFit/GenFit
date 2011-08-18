*
* $Id: gtau.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gtau.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE G3TAU
C.
C.    ******************************************************************
C.    *                                                                *
C.    *            Generates the process e+e- -> tau+ tau-             *
C.    *             can be called from GUKINE                          *
C.    *                                                                *
C.    *       Author    R.Hemingway  26 Sept 1985                      *
C.    *       Obsoleted in 3.15                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
      WRITE(CHMAIL,10000)
10000 FORMAT('Obsolete GTAU reached. This routine is now dummy.')
      CALL GMAIL(1,1)
C
      END
