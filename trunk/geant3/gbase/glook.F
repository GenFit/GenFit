*
* $Id: glook.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: glook.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GLOOK(NAME,IVECT,N,ILOOK)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Search position ILOOK of element NAME in array IVECT     *
C.    *       of length N                                              *
C.    *                                                                *
C.    *    ==>Called by : many GEANT and GEANG routines                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION IVECT(1)
      CHARACTER*(*) NAME
C.
C.    ------------------------------------------------------------------
C.
      CALL UCTOH(NAME,NAMEH,4,4)
      ILOOK=IUCOMP(NAMEH,IVECT,N)
C
      END
