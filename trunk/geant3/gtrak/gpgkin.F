*
* $Id: gpgkin.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gpgkin.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE G3PGKIN
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        Print kinematics of secondary tracks at the current     *
C.    *        interaction point.                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>,GUSTEP,G3DEBUG                        *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcking.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 UNITK,UNITE,UNITSK,UNITSE
C.
C.    ------------------------------------------------------------------
C.
      SEK=0.
      SET=0.
      DO 10 I=1,NGKINE
         ITY=GKIN(5,I)
         ET=GKIN(4,I)
         TOFDN=TOFD(I)*1.E+9
         JP=LQ(JPART-ITY)
         EK=ET-Q(JP+7)
         SEK=SEK+EK
         SET=SET+ET
         CALL G3EVKEV(EK,EKU,UNITK)
         CALL G3EVKEV(ET,ETU,UNITE)
         IF(I.LT.NGKINE)THEN
            WRITE(CHMAIL,1000)I,(IQ(JP+J),J=1,4),EKU,UNITK,ETU,UNITE
     +                         ,TOFDN
         ELSE
            CALL G3EVKEV(SEK,SEKU,UNITSK)
            CALL G3EVKEV(SET,SETU,UNITSE)
            WRITE(CHMAIL,2000)I,(IQ(JP+J),J=1,4),EKU,UNITK,ETU,UNITE
     +                         ,TOFDN,SEKU,UNITSK,SETU,UNITSE
         ENDIF
         CALL GMAIL(0,0)
  10  CONTINUE
*
 1000 FORMAT(4X,'===>',I3,1X,4A4,' EK=',F7.3,A,' ET=',F7.3,A,
     +       ' TOFD=',F10.3,' ns')
 2000 FORMAT(4X,'===>',I3,1X,4A4,' EK=',F7.3,A,' ET=',F7.3,A,
     +       ' TOFD=',F10.3,' ns   SEK=',F7.3,A,' SET=',F7.3,A)
      END
