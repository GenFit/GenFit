*
* $Id: gprsga.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gprsga.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:33  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE G3PRSGA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates cross-section in current material                  *
C.    *  for photon pair production and muon direct pair production.   *
C.    *                                                                *
C     *  semiempirical cross section formula of L.Urban is used        *
C     *  to estimate the photon mean free path in a given material     *
C     *            (see WRITEUP)                                       *
C.    *                                                                *
C.    *    ==>Called by : G3PHYSI                                      *
C.    *       Authors    R.Brun, G.Patrick, L.Urban  *********         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gccuts.inc"
C.
C.    ------------------------------------------------------------------
C.
C======>     A, Pair production by photons
C               ==========================
C               Cut on threshold kinetic energy. Special case for vacuum
C
      SST=0.
      IF(Z.LT.1.)GO TO 90
      EGAM   = ELOW(IEKBIN)
      IF (EGAM.LT.0.00099) GO TO 90
      IF (EGAM.LT.0.00199) EGAM=0.0014
C
      JPROB=LQ(JMA-4)
      IF(JMIXT.EQ.0)THEN
C
C             simple material (element)
C
          SST=G3PRSGG(Z,EGAM)/A
      ELSE
C
C             compound or mixture
C
         NLMAT=Q(JMA+11)
         NLM=IABS(NLMAT)
         DO 10 I=1,NLM
            II=JMIXT+NLM+I
            AA=Q(II-NLM)
            ZZ=Q(II)
            W=Q(II+NLM)
            SST=SST+W*G3PRSGG(ZZ,EGAM)/AA
  10     CONTINUE
      ENDIF
C
      SST=AVO*DENS*SST
C
  90  IF(SST.GT.0.)THEN
         Q(JPAIR+IEKBIN)=1./SST
      ELSE
         Q(JPAIR+IEKBIN)=BIG
      ENDIF
C
C======>     B, Direct pair production by muons
C               ===============================
C
      SST=0.
      IF(Z.LT.1.)GO TO 190
      T=ELOW(IEKBIN)
      IF(T.LT.1.)GO TO 190
      IF(T.LE.PPCUTM)GO TO 190
      IF(JMIXT.EQ.0)THEN
C
C            Element
C
         SST=G3PRSGM(Z,T,PPCUTM)
         SST=SST/A
      ELSE
C
C            Compound/Mixture
C
         NLMAT=Q(JMA+11)
         NLM  =IABS(NLMAT)
         SST  =0.
         DO 110 I=1,NLM
            II   = JMIXT+NLM+I
            AA   = Q(II-NLM)
            ZZ   = Q(II)
            WMAT = Q(II+NLM)
            SST  = SST+WMAT*G3PRSGM(ZZ,T,PPCUTM)/AA
 110     CONTINUE
      ENDIF
      SST=AVO*DENS*SST
C
 190  IF(SST.GT.0.)THEN
         Q(JPAIR+IEKBIN+NEK1)=1./SST
      ELSE
         Q(JPAIR+IEKBIN+NEK1)=BIG
      ENDIF
C
      END
