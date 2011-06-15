*
* $Id: randc.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: randc.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:43  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION FUNCTION RANDC(IDUMM)
#endif
#if defined(CERNLIB_SINGLE)
      FUNCTION RANDC(IDUMM)
#endif
C
      DIMENSION RND1(1)
C
      CALL GRNDM(RND1,1)
      RANDC = RND1(1)
      RETURN
      END
