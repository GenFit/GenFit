*
* $Id: gsdvs.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gsdvs.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:39  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE G3SDVS(CHNAME,CHMOTH,STEP,IAXIS,NUMED)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    To be replaced in the future by GSDVT                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
      CHARACTER*4 CHNAME,CHMOTH
      WRITE(CHMAIL,
     +'('' *** GSDVS: Obsolete Routine -- Call GSDVT instead'')')
      CALL GMAIL(1,0)
      WRITE(CHMAIL,
     +'(''            This code will be removed from next version'')')
      CALL GMAIL(0,1)
      CALL G3SDVT(CHNAME,CHMOTH,STEP,IAXIS,NUMED,0)
      END
