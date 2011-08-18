*
* $Id: gbrsga.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gbrsga.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:40  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE G3BRSGA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates cross-section in current material                  *
C.    *  for discrete(hard) electron BREMSSTRAHLUNG.                   *
C.    *                                                                *
C.    *    ==>Called by : G3PHYSI                                      *
C.    *       Authors    G.Patrick, L.Urban  *********                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcmulo.inc"
C.
C.    ------------------------------------------------------------------
C.
C=======>                 Electrons & Positrons
C
      SIG=0.
      IF(Z.LT.1.)GO TO 20
C
C             Calculate mean free path using total cross-section
C             formula of L.Urban (see Write-up)
C
      T      = ELOW(IEKBIN)
      IF(JMIXT.EQ.0)THEN
C
C           element
C
         SIG=G3BRSGE(Z,T,BCUTE)
         SIG=AVO*DENS*SIG/A
      ELSE
C
C           compound/mixture
C
         NLMAT=Q(JMA+11)
         NLM=IABS(NLMAT)
         SIG=0.
         DO 10 I=1,NLM
            J   = JMIXT+NLM+I
            AA  = Q(J-NLM)
            ZZ  = Q(J)
            WMAT= Q(J+NLM)
            S   = G3BRSGE(ZZ,T,BCUTE)
            S   = WMAT*S/AA
            SIG = SIG+AVO*DENS*S
  10     CONTINUE
      ENDIF
C
  20  IF(SIG.GT.0.)THEN
         Q(JBREM+IEKBIN)=1./SIG
         Q(JBREM+NEK1+IEKBIN)=1./(SIG*G3BFSIG(T,BCUTE))
      ELSE
         Q(JBREM+IEKBIN)=BIG
         Q(JBREM+NEK1+IEKBIN)=BIG
      ENDIF
C
C=======>                 Muons
C
      SIG=0.
      IF(Z.LT.1.)GO TO 40
C
C             Calculate mean free path using total cross-section
C             formula of L.Urban (see Write-up)
C
      IF(T.LT.1.)GO TO 40
      IF(JMIXT.EQ.0)THEN
C
C           element
C
         SIG=G3BRSGM(Z,T,BCUTM)
         SIG=AVO*DENS*SIG/A
      ELSE
C
C           compound/mixture
C
         NLMAT=Q(JMA+11)
         NLM=IABS(NLMAT)
         SIG=0.
         DO 30 I=1,NLM
            J   = JMIXT+NLM+I
            AA  = Q(J-NLM)
            ZZ  = Q(J)
            WMAT= Q(J+NLM)
            S   = G3BRSGM(ZZ,T,BCUTM)
            S   = WMAT*S/AA
            SIG = SIG+AVO*DENS*S
  30     CONTINUE
      ENDIF
C
  40  IF(SIG.GT.0.)THEN
         Q(JBREM+IEKBIN+2*NEK1)=1./SIG
      ELSE
         Q(JBREM+IEKBIN+2*NEK1)=BIG
      ENDIF
C
      END
