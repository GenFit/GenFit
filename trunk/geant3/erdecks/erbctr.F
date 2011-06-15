*
* $Id: erbctr.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: erbctr.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:35  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 15:37:34  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERBCTR (X)
*
C.    ******************************************************************
C.    *                                                                *
C.    *       After the back-tracking transforms                       *
C.    *                  the transport matrix back to the original     *
C.    *                  particle direction                            *
C.    *                                                                *
C.    *       Arguments                                                *
C.    *           X        - Transport Matrix                          *
C.    *                                                                *
C.    *    ==>Called by : ERPEND                                       *
C.    *       Author    E.Nagy  *********                              *
C.    *                                                                *
C.    ******************************************************************
*
         DOUBLE PRECISION X(5,5)
*
      DO 10 I = 1, 5
         DO 5 K = 1, 5
            IF (I.EQ.2) X(I,K) = -X(I,K)
            IF (K.EQ.2) X(I,K) = -X(I,K)
            IF (I.EQ.4) X(I,K) = -X(I,K)
            IF (K.EQ.4) X(I,K) = -X(I,K)
    5    CONTINUE
   10 CONTINUE
*                                                           END ERBCTR
      END
