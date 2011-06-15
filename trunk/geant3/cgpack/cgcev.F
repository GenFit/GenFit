*
* $Id: cgcev.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cgcev.F,v $
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
      SUBROUTINE CGCEV(IFVIS,CG)
************************************************************************
*                                                                      *
*     Name: CGCEV                                                      *
*     Author: E. Chernyaev                       Date:    05.02.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Change edge visibility                                 *
*                                                                      *
*     References: CGTSTR                                               *
*                                                                      *
*     Input:   IFVIS - flag (1 - visible, 0 - unvisible, -1 - inverse) *
*     Output: CG(*)  - CG-object                                       *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      CG(*)
*-
      CALL CGTSTR(CG,IREP)
      IF (IREP .LT. 0)          GOTO 999
      NFACE  = CG(KCGNF)
      IF (NFACE .EQ. 0)         GOTO 999
      JCG    = LCGHEA
      DO 200 NF=1,NFACE
        NEDGE  = CG(JCG+KCGNE)
        JCG    = JCG + LCGFAC
        DO 100 NE=1,NEDGE
          IF (IFVIS .GT. 0)     CG(JCG+KCGAE) = 0.
          IF (IFVIS .EQ. 0)     CG(JCG+KCGAE) =-1.
          IF (IFVIS .LT. 0)     CG(JCG+KCGAE) =-(CG(JCG+KCGAE) + 1.)
          JCG    = JCG + LCGEDG
  100     CONTINUE
  200   CONTINUE
*
  999 RETURN
      END
