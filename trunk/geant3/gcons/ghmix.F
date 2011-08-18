*
* $Id: ghmix.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ghmix.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE G3HMIX(A, W, N, AEFF)
C
C     ******************************************************************
C     *                                                                *
C     *   Works out an effective atomic weight AEFF for a material     *
C     * with N elements of atomic weight A(I) in proportion W(I) by    *
C     * weight.  The criterion is that the hadronic interaction        *
C     * length of a 5 GeV/c pion is correct.  Errors on the calculated *
C     * hadronic interaction length for other momenta and other        *
C     * particles in GEANT version 3.04 are less than 1% in most cases.*
C     * For details see Memorandum OPAL/0037N/JA/md, ref. Hadron       *
C     * Milestone 84/003, Calculation of Hadronic Interaction Lengths  *
C     * for Mixtures.                                                  *
C.    *                                                                *
C.    *    ==>Called by : G3SMIXT                                      *
C.    *       Author    J.Allison  *********                           *
C     *                                                                *
C     ******************************************************************
C
#include "geant321/gcunit.inc"
      DIMENSION A(N), W(N)
C.
C.    ------------------------------------------------------------------
C.
C         Work out PINT which is proportional to the interaction
C         probability.  Also work out the mean atomic weight, i.e. that
C         weighted by proportion by numbers, as starting value for
C         iterative method of finding AEFF.
C
      PINT   = 0.
      AINV   = 0.
      WTOT   = 0.
      DO 10 I = 1, N
         PINT   = PINT + W(I) * G3HSIGM(5., 8, A(I)) / A(I)
         AINV   = AINV + W(I) / A(I)
         WTOT   = WTOT + W(I)
  10  CONTINUE
C
      IF ( ABS ( WTOT - 1. ) .GT. 0.01 ) GO TO 98
C
C         Work out AEFF which gives PINT for 5 GeV/c pion.
C         (This is a short Newton's method loop.)
C
      AEFF   = 1. / AINV
      PNEW   = G3HSIGM ( 5., 8, AEFF ) / AEFF
      DA     = 1.
      ITER   = 0
   20 CONTINUE
      ITER   = ITER + 1
      AEFF   = AEFF + DA
      POLD   = PNEW
      DAOLD  = DA
      PNEW   = G3HSIGM ( 5., 8, AEFF ) / AEFF
      DP     = PNEW - POLD
      DA     = (PINT - PNEW ) * DAOLD / DP
      IF ( ( ITER .GT. 1 ) .AND. ( ABS ( DA ) .GT. ABS ( DAOLD ) ) )
     +                            GO TO 97
      IF ( ABS ( DA ) .GT. 0.01 ) GO TO 20
C
      RETURN
C
C            Error conditions.
C
  97  CONTINUE
      WRITE (CHMAIL,197)
      CALL GMAIL(0,0)
 197  FORMAT ( ' ***** GHMIX : ',
     +'HADRONIC INTERACTION MIXTURE ROUTINE NOT CONVERGING')
      RETURN
C
  98  CONTINUE
      WRITE (CHMAIL,198) WTOT
      CALL GMAIL(0,0)
 198  FORMAT ( ' ***** GHMIX : ',
     +'FRACTIONS BY WEIGHT OF MIXTURES DO NOT ADD UP TO 1',F10.4)
      END
