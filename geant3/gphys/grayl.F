*
* $Id: grayl.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: grayl.F,v $
* Revision 1.2  2003/11/28 11:23:56  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.2  1996/11/06 12:53:21  ravndal
* Correct the rotation into the GEANT system
*
* Revision 1.1.1.1  1995/10/24 10:21:33  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE G3RAYL
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Generates COHERENT SCATTERING mechanism                       *
C.    *                                                                *
C.    *    ==>Called by : G3TGAMA                                      *
C.    *       Author    G.Tromba, P.Bregant                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcking.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
      DIMENSION VROT(3)
      DIMENSION RNDM(2)
      LOGICAL ROTATE
C.
C.    ------------------------------------------------------------------
C.
      KCASE = NAMEC(25)
      NGKINE=0
*
      IMAX=NEKBIN
      DO 2 I=2,NEK1
         IF (Q(JRAYL+I).GT.0.9*BIG) THEN
            IMAX = I-1
            GO TO 2
         ENDIF
   2  CONTINUE
C
C             Extract a value for the momentum transfer according
C             with the FORM FACTORS distribution.
C
      IF (IEKBIN.LT.IMAX) THEN
        GEKRT1=1.-GEKRAT
        A2MAX=GEKRT1*Q(JRAYL+NEK1+IEKBIN)+GEKRAT*Q(JRAYL+NEK1+IEKBIN+1)
      ELSE
        A2MAX=Q(JRAYL+NEK1+IEKBIN)
      ENDIF
      GEKIN1=1./GEKIN
   5  CALL GRNDM(RNDM,2)
      A2=RNDM(1)*A2MAX
      DO 10 I=1,(IMAX-1)
         IF (A2.LE.Q(JRAYL+I+NEK1)) THEN
            IE=I
            GO TO 20
         ENDIF
10    CONTINUE
      IE=IMAX
20    IF (IE.GT.1) THEN
        GA2RAT =(A2-Q(JRAYL+NEK1+IE-1))/(Q(JRAYL+NEK1+IE)
     +          -Q(JRAYL+NEK1+IE-1))
        QNEW   = ELOW(IE-1)*(1.-GA2RAT)+ELOW(IE)*GA2RAT
      ELSE
        GA2RAT  = A2/Q(JRAYL+NEK1+IE)
        QNEW    = ELOW(IE)*GA2RAT
      ENDIF
C
C            Find out cos(theta)
      COSTH=1.-0.5*(QNEW*GEKIN1)**2
      IF (COSTH.LT.-1.) COSTH=-1.
      IF (COSTH.GT.1.)  COSTH=1.
      IF(RNDM(2).GT.0.5*(1.+COSTH**2)) GO TO 5
      SINTH=SQRT((1.-COSTH)*(1.+COSTH))
C
C            Extract PH  (PHI unif. between 0 and 2PI)
      CALL GRNDM(RNDM,2)
      PH=RNDM(1)*TWOPI
      COSPH=COS(PH)
      SINPH=SIN(PH)
C
C            Set new VECT(I), i=4,5,6
      VROT(1)=GEKIN*SINTH*COSPH
      VROT(2)=GEKIN*SINTH*SINPH
      VROT(3)=GEKIN*COSTH
C
C             Rotate scattered photon into GEANT system
C
      CALL G3FANG(VECT(4),COSTH,SINTH,COSPH,SINPH,ROTATE)
      IF (ROTATE) CALL G3DROT(VROT,COSTH,SINTH,COSPH,SINPH)
 
      VMROT=1./SQRT(VROT(1)*VROT(1)+VROT(2)*VROT(2)+VROT(3)*VROT(3))
      VECT(4)=VROT(1)*VMROT
      VECT(5)=VROT(2)*VMROT
      VECT(6)=VROT(3)*VMROT
C
C            Update probabilities
C
*
      ZINTRA=-LOG(RNDM(2))
      SLRAYL=SLENG
      STEPRA=BIG
*
C
      END
 
