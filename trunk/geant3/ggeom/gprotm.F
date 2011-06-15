*
* $Id: gprotm.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gprotm.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:39  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE G3PROTM(NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print rotation matrices structure JROTM       *
C.    *       NUMB     Rotation matrix number                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *         Author  R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF (JROTM.LE.0)  GO TO 99
C
      IF (NUMB  .EQ.0)THEN
         WRITE (CHMAIL,1000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = IQ(JROTM-2)
      ELSE
         N1     = ABS(NUMB)
         N2     = ABS(NUMB)
      ENDIF
C
      IF(NUMB.GE.0) THEN
         WRITE (CHMAIL,1001)
         CALL GMAIL(0,1)
      ENDIF
C
      DO 30 I=N1,N2
      JR = LQ(JROTM-I)
      IF (JR.LE.0)   GO TO 30
C
      IFLAG = Q(JR+10)
      WRITE (CHMAIL,1002) I,(Q(JR+J),J = 1,9),IFLAG,(Q(JR+J),J = 11,16)
      CALL GMAIL(0,0)
   30 CONTINUE
C
 1000 FORMAT ('0',51('='),1X,'ROTATION MATRICES',2X,50('='))
 1001 FORMAT ('0','MATRIX',1X,'<',17('-'),' ROTATION MATRIX',1X,17('-')
     +,'>',' FLAG','  THET1','   PHI1','  THET2','   PHI2','  THET3'
     +,'   PHI3')
 1002 FORMAT (' ',I6,9F6.3,I5,6F7.2)
   99 RETURN
      END
