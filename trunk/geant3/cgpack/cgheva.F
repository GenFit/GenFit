*
* $Id: cgheva.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cgheva.F,v $
* Revision 1.1.1.1  2002/06/16 15:17:54  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:04  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGHEVA(H,HISI)
*SG
#include "geant321/cghpar.inc"
*
      REAL H(*)
*-
      HISI=0
      NHSIZE = H(KHSIZE)
*      NT     = H(KHNT)
      NXYZ   = H(KHNXYZ)
      NFACE  = H(KHNFAC)
      NEDGE  = H(KHNEDG)
*           S P A C E   F O R   N O D E S
      NX     = 3*NXYZ
*           S P A C E   F O R   F A C E S
      NF     = (NFACE+NEDGE*2)/I4SIZE + 1
      NP     = NFACE
      NDF    = 6*NFACE
*           S P A C E   F O R   T R E E
      NTR    = 4*NFACE/I4SIZE
      NST    = NFACE/I4SIZE + 1
      NAL    = NFACE
      NAR    = NFACE
*            C H E C K   S P A C E
      NTOTAL = LHHEAD + NX + NF + NP + NDF + NTR + NST + NAL + NAR
      IF(NTOTAL .GT. NHSIZE)THEN
        HISI=NTOTAL+1
      ENDIF
 999  END
*SG
