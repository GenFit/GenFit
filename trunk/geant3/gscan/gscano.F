*
* $Id: gscano.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gscano.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:36  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GSCANO
*
************************************************************************
*             SCAN utility called at the end of each event.            *
*                                                                      *
*          ==>Called by : <USER>,GUKINE                                *
*             Author    R.Brun  *********                              *
************************************************************************
*                                                                      *
#include "geant321/gcbank.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcscan.inc"
C.
C.    ------------------------------------------------------------------
C.
*
*          Insert the temporary LSTEMP structure into LSCAN
*
      IF(IPHI.GT.IQ(LSCAN-2))GO TO 10
      IF(LQ(LSCAN-IPHI).NE.0)GO TO 10
      NPUSH=ISCUR-IQ(LSTEMP-1)-1
      IQ(LSTEMP-5)=IPHI
      CALL MZPUSH(IXCONS,LSTEMP,0,NPUSH,'I')
      CALL ZSHUNT(IXCONS,LSTEMP,LSCAN,-IPHI,0)
      GO TO 999
   10 IF(LSTEMP.NE.0)CALL MZDROP(IXCONS,LSTEMP,' ')
      LSTEMP=0
*
  999 END
