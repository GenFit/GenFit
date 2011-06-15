*
* $Id: gannii.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gannii.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE G3ANNII
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates cross-section in current material                  *
C.    *  for positron annihilation                                     *
C.    *                                                                *
C.    *    ==>Called by : G3PHYSI                                      *
C.    *       Authors    R.Brun, L.Urban  *********                    *
C.    *       10/06/93: modified by Georges Azuelos (Vancouver)        *
C     *                 to include 1-quantum annihilation              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
C.
C.    ------------------------------------------------------------------
C.
      SIG=0.
      IF(Z.GE.1.) THEN
         XE=ELOW(IEKBIN)+EMASS
         GAM=XE/EMASS
         GAM2=GAM**2
         GAM1=MAX(GAM2-1.,0.)
         C=SQRT(GAM1)
*
         SIG=(GAM2+4.*GAM+1.)*LOG(GAM+C)/GAM1-(GAM+3.)/C
         SIG=0.5*Q(JPROB+17)*SIG/(GAM+1.)
         PE=SQRT((XE-EMASS)*(XE+EMASS))
*
         GVE=PE/EMASS
         SIG1=GAM2+2.*(GAM+2.)/3.-(GAM+2.)/GVE*LOG(GAM+GVE)
         SIG1=2.*Q(JPROB+18)*SIG1/(GVE*(GAM+1.)**2)
*
         SIG=SIG+SIG1
      ENDIF
      IF(SIG.GT.0.)THEN
         Q(JANNI+IEKBIN)=1./SIG
      ELSE
         Q(JANNI+IEKBIN)=BIG
      ENDIF
C
      END
