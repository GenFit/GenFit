*
* $Id: gstren.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gstren.F,v $
* Revision 1.2  2003/11/28 11:23:56  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      FUNCTION G3STREN (GAMMA,ECUT,STEP)
C.    ******************************************************************
C.    *                                                                *
C.    *   Simulation of energy loss straggling in thin layers.         *
C.    *   Sampling is done from the tables which are prepared          *
C.    *   in GSTINI.                                                   *
C.    *                                                                *
C.    *   ==> Called by : G3TELEC, G3THADR, G3TMUON                    *
C.    *       Authors   : K. Lassila-Perini, I. Gavrilenko             *
C.    *                                                                *
C.    ******************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcstra.inc"
C
      DIMENSION ANGL(200),RNDM(2)
C...
*   Reading the tables, if the medium has not changed
      NP     = Q(JTSTRA+1)
*
      EC       = 1.E9*ECUT
      GSTLOC   = 0.
      IF(STEP.LE.0.) GO TO 90
C...
*   Interpolate the tables for the current gamma-factor
      GAML     = LOG(GAMMA)
      IF(GAML.GT.GAMLOG(21)) THEN
         DO 10 J=1,NP
            ANGL(J)=Q(JTSTCO+J*21)
   10    CONTINUE
      ELSE
C...
         DO 20 I = 2,21
            IF(GAML.LE.GAMLOG(I)) THEN
               DGG = (GAML-GAMLOG(I-1))/(GAMLOG(I)-GAMLOG(I-1))
               IP=I
               GO TO 30
            ENDIF
   20    CONTINUE
   30    DG1=1.-DGG
C...
         DO 40 K = 1,NP
            ANGL(K) = Q(JTSTCO+(K-1)*21+IP-1)*DG1+
     +                Q(JTSTCO+(K-1)*21+IP  )*DGG
   40    CONTINUE
      ENDIF
*
*    The sampling from the tables
      ANMIN   = EXP(ANGL(NP))
      ANMAX   = EXP(ANGL( 1))
*
      DAN      = ANMAX-ANMIN
      AMU      = DAN*STEP
      CALL G3POISS(AMU,N,1)
*      NICOLL   = MAX(N,1)
      NICOLL = N
      DO 80  I  = 1,NICOLL
   50    K1 = 1
         K2 = NP
         CALL GRNDM(RNDM,1)
         R = LOG(RNDM(1)*DAN+ANMIN)
   60    K = (K1+K2)/2
         IF(R.GT.ANGL(K)) THEN
            K2=K
         ELSEIF(R.LT.ANGL(K)) THEN
            K1=K
         ELSE
            E=Q(JTSTEN+K)
            GOTO 70
         ENDIF
         IF(K2-K1.GT.1 ) GO TO 60
*
         E = Q(JTSTEN+K2)+(R-ANGL(K2))*
     +       (Q(JTSTEN+K1)-Q(JTSTEN+K2))/(ANGL(K1)-ANGL(K2))
   70    EE = EXP(E)
         IF(EE.GT.EC) GO TO 50
*
* *** Total energy calculation
         GSTLOC = GSTLOC+EE
   80 CONTINUE
*
   90 G3STREN  = GSTLOC*1.E-9
      END
