*
* $Id: ginpar.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ginpar.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:39  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINPAR(X,PAR,IYES)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Routine to determine whether a point X is inside the        *
C.    *    parallelepiped with parameters PAR. IYES is returned as 1   *
C.    *    if it is; as zero if it isn't.                              *
C.    *                                                                *
C.    *    ==>Called by : GINME                                        *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X(3),PAR(6)
C.
C.           -------------------------------------
C.
      IYES=0
      IF(ABS(X(3)).GT.PAR(3)) GO TO 99
C
      YT=X(2)-PAR(6)*X(3)
      IF(ABS(YT).GT.PAR(2)) GO TO 99
C
      XT=X(1)-PAR(5)*X(3)-PAR(4)*YT
      IF(ABS(XT).GT.PAR(1)) GO TO 99
C
      IYES=1
C
   99 CONTINUE
      END
