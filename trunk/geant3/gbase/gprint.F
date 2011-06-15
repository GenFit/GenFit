*
* $Id: gprint.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gprint.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3PRINT(CHNAME,NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print data structures                         *
C.    *                                                                *
C.    *       CHNAME   name of a data structure                        *
C.    *       NUMB     data structure number                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT  ,<GXINT> GINC4               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      CHARACTER*4 KNAMES(11),NAME
      DIMENSION NHAMES(11)
      CHARACTER*(*) CHNAME
      SAVE IFIRST,KNAMES,NHAMES
      DATA IFIRST/0/
      DATA KNAMES/'DIGI','JXYZ','HITS','KINE','MATE','VOLU'
     +           ,'ROTM','SETS','TMED','PART','VERT'/
C.
C.    ------------------------------------------------------------------
C.
      NAME=CHNAME
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL UCTOH(KNAMES,NHAMES,4,44)
      ENDIF
C
      CALL GLOOK(NAME,NHAMES,11,JUMP)
C
      IF(JUMP.EQ. 1) CALL G3PDIGI ('*','*')
      IF(JUMP.EQ. 2) CALL G3PJXYZ (NUMB)
      IF(JUMP.EQ. 3) CALL G3PHITS ('*','*')
      IF(JUMP.EQ. 4) CALL G3PKINE (NUMB)
      IF(JUMP.EQ. 5) CALL G3PMATE (NUMB)
      IF(JUMP.EQ. 6) CALL G3PVOLU (NUMB)
      IF(JUMP.EQ. 7) CALL G3PROTM (NUMB)
      IF(JUMP.EQ. 8) CALL G3PSETS ('*','*')
      IF(JUMP.EQ. 9) CALL G3PTMED (NUMB)
      IF(JUMP.EQ.10) CALL G3PPART (NUMB)
      IF(JUMP.EQ.11) CALL G3PVERT (NUMB)
C
      END
