*
* $Id: granor.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: granor.F,v $
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
      SUBROUTINE GRANOR(A,B)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       To generate 2 numbers A and B following a NORMAL         *
C.    *       distribution (mean=0 sigma=1.)                           *
C.    *         Copy of the CERN Library routine RANNOR                *
C.    *                                                                *
C.    *    ==>Called by : <USER>, many GEANT routines                  *
C.    *       Author    F.Carminati *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION RNDM(2)
*
      CALL GRNDM(RNDM,2)
      Y=RNDM(1)
      Z=RNDM(2)
      X=6.283185*Z
      A1=SQRT (-2.0*LOG(Y))
      A=A1*SIN (X)
      B=A1*COS (X)
      RETURN
      END
 
