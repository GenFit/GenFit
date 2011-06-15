*
* $Id: ghclos.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: ghclos.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3HCLOS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        SUBR. G3HCLOS                                           *
C.    *                                                                *
C.    *   For every sensitive detector, calls GGDETV in order to       *
C.    *   complete the JD bank with the list of volume names which     *
C.    *   permit to identify uniquely a given physical detector, the   *
C.    *   list of bit numbers to pack the corresponding volume copy    *
C.    *   numbers,  and the generic path(s) through the JVOLUM tree.   *
C.    *                                                                *
C.    *   Called by : G3GCLOS, <USER>                                  *
C.    *   Author    : F.Bruyant                                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
C.
C.    ------------------------------------------------------------------
*
      IF (JSET.NE.0) THEN
*
* ***   Scan the volume structure to retrieve the path through
*       the physical tree for all sensitive detectors
*
        NSET = IQ(JSET-1)
        DO 20 ISET=1,NSET
          JS = LQ(JSET-ISET)
          NDET = IQ(JS-1)
          DO 10 IDET=1,NDET
            CALL GGDETV (ISET, IDET)
   10     CONTINUE
   20   CONTINUE
      ENDIF
*                                                             END GHCLOS
      END
