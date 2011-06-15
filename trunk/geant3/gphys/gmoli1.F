*
* $Id: gmoli1.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gmoli1.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:25  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.2  2002/05/13 12:40:58  hristov
* Dummy subroutines to avoid files with no code in
*
* Revision 1.1.1.1  1999/05/18 15:55:20  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_VER314)
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE G3MOLI1(AC,ZC,DENS,BETA2,OMC)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Initializes material constants OMC,CHC for Moliere         *
C.    *     scattering (case of a simple material)                     *
C.    *                                                                *
C.    *    ==>Called by : G3MULOF                                      *
C.    *         Author M.S. Dixit NRCC Ottawa    *********             *
C.    *                                                                *
C.    ******************************************************************
C.
C     ZETA=factor for scattering by atomic electrons
C
      DATA ZETA / 1. /
*
*     ------------------------------------------------------------------
*
      ZS=(1./AC)*ZC*(ZC+ZETA)
      ZE=-0.66666*ZS*LOG(ZC)
      ZX=ZS*LOG(1.+3.34*(ZC/(BETA2*137.036))**2)
      OMC=6702.33*ZS*DENS*EXP((ZE-ZX)/ZS)
*
      END
#else
      SUBROUTINE GMOLI1_DUMMY
      END
#endif
