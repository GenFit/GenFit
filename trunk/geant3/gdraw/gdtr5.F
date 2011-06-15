*
* $Id: gdtr5.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gdtr5.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:03  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR5(JM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Build one level of optimized tree                        *
C.    *                                                                *
C.    *       JM = mother node (input)                                 *
C.    *                                                                *
C.    *    ==>Called by : GDTR6                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      JM1=IQ(JSCA+JM)
      JX=IQ(JXON+JM)
      JX1=IQ(JSCA+JX)
C
C             Copy multiplicity from non-optimized to optimized tree
C
      IF (IQ(JDVM+JM).GT.0) THEN
         IQ(JDUP1+JM1)=-IQ(JDVM+JM)
      ELSE
***     IQ(JDUP1+JM1)=IQ(JDUP+JM)*IQ(JPSM+JM)
         IQ(JDUP1+JM1)=IQ(JDUP+JM)
      ENDIF
C
      IF (JX.EQ.0) THEN
         IQ(JXON1+JM1)=0
      ELSE
         IQ(JXON1+JM1)=JX1
         LXON=IQ(JNAM+JX)
         IQ(JNAM1+JX1)=LXON
         IQ(JMOT1+JX1)=JM1
         LBRO=IQ(JBRO+JX)
   10    IF (LBRO.EQ.0) GO TO 20
         IF (IQ(JDUP+LBRO).NE.0) THEN
            LBRO1=IQ(JSCA+LBRO)
            LXON=IQ(JNAM+LBRO)
            IQ(JBRO1+JX1)=LBRO1
            IQ(JNAM1+LBRO1)=LXON
            IQ(JMOT1+LBRO1)=JM1
            JX1=LBRO1
         ENDIF
         LBRO=IQ(JBRO+LBRO)
         GO TO 10
   20    CONTINUE
         IQ(JBRO1+JX1)=0
      ENDIF
      END
