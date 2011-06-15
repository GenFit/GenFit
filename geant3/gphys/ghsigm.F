*
* $Id: ghsigm.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ghsigm.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:36  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      FUNCTION G3HSIGM(P,IPART,A)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Returns absorption cross-section in millibarns               *
C.    *   for a particle with                                          *
C.    *                momentum P (GEV/C)                              *
C.    *                GEANT type IPART                                *
C.    *      on a nucleus of atomic number A                           *
C.    *                                                                *
C.    *    ==>Called by : G3PHADR                                      *
C.    *       Author    R.Barlow  *********                            *
C.    *       Modified to deal with K0 and extend to 41 particle       *
C.    *                types by John Allison, 31/1/84. *               *
C.    *                                                                *
C.    ******************************************************************
C.
C              The array NTYPE maps the GEANT particle type (IPART)
C              onto G3HSIG's internal particle type.  Simple assumptions
C              are made for highly flavoured particles.  The K0 is
C              given a non-existent internal type 6 which causes the
C              mean of the K+ and K- cross sections to be calculated.
C
      DIMENSION NTYPE(41)
      DATA NTYPE/6*0,3*5,6,3,4,2*1,2,6,
     +  8*1,8*2,2*0,4,3,2*6,2*3,1/
C.
C.    ------------------------------------------------------------------
C.
      G3HSIGM = 1.E-20
      IF (IPART .GT. 41) RETURN
      ITYPE  = NTYPE(IPART)
      IF (ITYPE .EQ. 0) RETURN
C
      IF (ITYPE .NE. 6)THEN
         G3HSIGM = G3HSIG(P, ITYPE, A)
      ELSE
         G3HSIGM = 0.5 * (G3HSIG(P, 3, A) + G3HSIG(P, 4, A))
      ENDIF
C
      END
