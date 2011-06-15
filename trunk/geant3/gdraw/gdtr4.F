*
* $Id: gdtr4.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gdtr4.F,v $
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
      SUBROUTINE GDTR4
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Collect multiplicity for all JVOLUM structure            *
C.    *                                                                *
C.    *    ==>Called by : GDTR0                                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      DO 10 J=1,NUMNOD
         IQ(JSCA+J)=0
         IQ(JDUP+J)=0
   10 CONTINUE
      IQ(JDUP+1)=1
      DO 30 J=1,NUMNOD
         IF (IQ(JSCA+J).EQ.0) THEN
            CALL GDTR3(J)
            DO 20 I=J,NUMNOD
               IF (IQ(JSCA+I).EQ.-1) IQ(JSCA+I)=0
   20       CONTINUE
         ENDIF
   30 CONTINUE
      END
