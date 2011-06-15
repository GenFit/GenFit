*
* $Id: guview.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: guview.F,v $
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
      SUBROUTINE GUVIEW(IDS,IVS,ICS,IVIEW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine for interactive version                     *
C.    *                                                                *
C.    *    ==>Called by : <GXINT>, GINC1                               *
C.    *                                                                *
C.    ******************************************************************
C.
C.
C.    ------------------------------------------------------------------
C.
      CHARACTER*(*) ICS
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#if !defined(CERNLIB_USRJMP)
      CALL GUSTEP
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT0(JUSTEP)
#endif
      END
