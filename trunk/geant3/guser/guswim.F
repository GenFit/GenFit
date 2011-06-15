*
* $Id: guswim.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: guswim.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUSWIM(CHARGE,STEP,VECT,VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to control tracking of one track            *
C.    *       in a magnetic field                                      *
C.    *                                                                *
C.    *    ==>Called by : G3TELEC,G3THADR,G3TMUON                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gctmed.inc"
      DIMENSION VECT(6),VOUT(6)
C.
C.    ------------------------------------------------------------------
C.
         IF (IFIELD.EQ.3)THEN
            CALL G3HELX3(FIELDM*CHARGE,STEP,VECT,VOUT)
         ELSEIF(IFIELD.EQ.2)THEN
            CALL G3HELIX(CHARGE,STEP,VECT,VOUT)
         ELSE
            CALL G3RKUTA(CHARGE,STEP,VECT,VOUT)
         ENDIF
C
      END
