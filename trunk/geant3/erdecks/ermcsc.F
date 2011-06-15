*
* $Id: ermcsc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ermcsc.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:35  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERMCSC (E)
*
C.    ******************************************************************
C.    *                                                                *
C.    *       Calculates Multiple Scattering Matrix in a step          *
C.    *                                                                *
C.    *           E        - Error Matrix (Triangular)                 *
C.    *           STEP     - Step length in cm                         *
C.    *           RADL     - Radiation length in cm                    *
C.    *           VECT(7)  - Momentum in GeV                           *
C.    *           CLA      - Cosine of the deep angle                  *
C.    *                                                                *
C.    *    ==>Called by : ERPROP                                       *
C.    *       Author    E.Nagy  *********                              *
C.    *                                                                *
C.    ******************************************************************
*
#include "geant321/gcmate.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
*
      DOUBLE PRECISION E(15)
      REAL*4 RADLOCAL
*
      PBETA = VECT(7)*VECT(7)/GETOT
      
      CRLOG = 1.
CC      CRLOG = 1. + 0.201*LOG(RI)
CC      IF (CRLOG.LE.0.) CRLOG = 0.
CC      D     = 1.8496E-4*RI*((CHARGE*CRLOG/PBETA)**2)
*
* modified by Andrea Fontana and Alberto Rotondi to include 
* the correct variance for Coulomb multiple scattering
* (see Fruhwirth and Regler, NIM A 456 (2001) 369)
* march 2007
*
* begin of change
*
      RADLOCAL = (Z+1)*LOG(287.*Z**(-0.5))/(Z*LOG(159.*Z**(-0.333333)))
      RI    = STEP/(RADL*RADLOCAL)
      D     = 2.25E-4*RI*((CHARGE*CRLOG/PBETA)**2)
*
* end of change
*
      S1    = D*STEP**2/3.
      S2    = D
      S3    = D*STEP/2.
      CLA = SQRT (VECT(4)**2+VECT(5)**2)
      E(6)  = E(6)+S2
      E(9)  = E(9)-S3
      E(10) = E(10)+S2/CLA**2
      E(11) = E(11)+S3/CLA
      E(13) = E(13)+S1
      E(15) = E(15)+S1
*
      END
