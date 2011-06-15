*
* $Id: gukine.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gukine.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE GUKINE
*
************************************************************************
*                                                                      *
*             Generates Kinematics for primary tracks                  *
*                                                                      *
************************************************************************
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcomis.inc"
*
      DIMENSION VERTEX(3),PLAB(3)
      DIMENSION RNDM(2)
*
*     -----------------------------------------------------------------
*
#if defined(CERNLIB_COMIS)
      IF(JUKINE.NE.0)THEN
         CALL CSJCAL(JUKINE,0,X,X,X,X,X,X,X,X,X,X)
         GO TO 99
      ENDIF
#endif
 
      VERTEX(1)=VSCAN(1)
      VERTEX(2)=VSCAN(2)
      VERTEX(3)=VSCAN(3)
      IF(IKINE.GT.100)THEN
         IK=IKINE-100
         THETA=PKINE(2)*DEGRAD
         PHI=PKINE(3)*DEGRAD
      ELSE
         IK=IKINE
         CALL GRNDM(RNDM,2)
         THETA=PI*RNDM(1)
         PHI=TWOPI*RNDM(2)
      ENDIF
C
      PLAB(1) = PKINE(1)*SIN(THETA)*COS(PHI)
      PLAB(2) = PKINE(1)*SIN(THETA)*SIN(PHI)
      PLAB(3) = PKINE(1)*COS(THETA)
C
      CALL G3SVERT(VERTEX,0,0,0,0,NVERT)
      CALL G3SKINE(PLAB,IK,NVERT,0,0,NT)
*
*              Kinematic debug (controled by ISWIT(1))
*
      IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
        CALL G3PRINT('VERT',0)
        CALL G3PRINT('KINE',0)
      ENDIF
*
  99  END
