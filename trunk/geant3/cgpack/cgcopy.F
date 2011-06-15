*
* $Id: cgcopy.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cgcopy.F,v $
* Revision 1.1.1.1  2002/06/16 15:17:54  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:04  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGCOPY(A,NMAX,B)
************************************************************************
*                                                                      *
*     Name: CGCOPY                                                     *
*     Author: E. Chernyaev                       Date:    02.11.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Copy one CG-object to another                          *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: A(*)  - 1-st CG-object                                    *
*            NMAX  - max length B array                                *
*                                                                      *
*     Output: B(*) - resulting CG-object                               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      A(*),B(*)
*-
      NSIZE     = A(KCGSIZ)
      B(KCGSIZ) = 0.
      IF (NSIZE .LE. 0)         GOTO 999
      IF (NMAX .LE. 0)          GOTO 999
      IF (NSIZE .GT. NMAX)      B(KCGSIZ) = NMAX - NSIZE
      IF (NSIZE .GT. NMAX)      GOTO 999
      DO 100 I=1,NSIZE
        B(I) = A(I)
  100   CONTINUE
*
  999 RETURN
      END
