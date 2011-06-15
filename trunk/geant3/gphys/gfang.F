*
* $Id: gfang.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfang.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 06/07/94  18.28.55  by  S.Giani
*-- Author :
      SUBROUTINE G3FANG(P,COSTH,SINTH,COSPH,SINPH,ROTATE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Finds theta-phi angles for a particle with director cosines P *
C.    *  The logical flag ROTATE indicates whether the the director    *
C.    *  cosines are along the Z axis (.FALSE.) or not (.TRUE.).       *
C.    *                                                                *
C.    *    ==>Called by :                                              *
C.    *       Author    F.Carminati  *********                         *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION P(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION DUX, DUY, DUZ, DSITH2, ONE, DNORM
#endif
      LOGICAL ROTATE
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      DUX      = P(1)
      DUY      = P(2)
      DUZ      = P(3)
      ROTATE   = .TRUE.
      IF(ABS(DUZ).GE.0.85) THEN
         DSITH2=DUX**2+DUY**2
         IF(DSITH2.GT.0.) THEN
            COSTH  = SIGN(ONE,DUZ)*SQRT(ONE-DSITH2)
            DSITH  = SQRT(DSITH2)
            SINTH  = DSITH
            COSPH  = DUX/DSITH
            SINPH  = DUY/DSITH
         ELSEIF(DUZ.GT.0.)THEN
            ROTATE = .FALSE.
            COSTH =  1.
            SINTH =  0.
            COSPH =  1.
            SINPH =  0.
         ELSE
            COSTH = -1.
            SINTH =  0.
            COSPH =  1.
            SINPH =  0.
         ENDIF
      ELSE
         COSTH  = DUZ
         DSITH  = SQRT((ONE+DUZ)*(ONE-DUZ))
         SINTH  = DSITH
         DNORM  = ONE/SQRT(DUX**2+DUY**2)
         COSPH  = DUX*DNORM
         SINPH  = DUY*DNORM
      ENDIF
*
      END
