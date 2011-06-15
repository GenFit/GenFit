*
* $Id: rfinv.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: rfinv.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:48  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:33  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:03  mclareni
* Kernlib
*
*
#include "kernnum/pilot.h"
          SUBROUTINE          RFINV(N,A,IDIM,IR)
          INTEGER             IR(*)
          REAL                A(IDIM,*), X, Y, TI
          CHARACTER*6         HNAME
#if !defined(CERNLIB_NUMRDBLE)
          REAL                ZERO, S31, S32, S33, S34, DOTF
          DOTF(X,Y,S31)  =  X*Y + S31
          DATA      ZERO      /  0.  /
#endif
#if defined(CERNLIB_NUMRDBLE)
          DOUBLE PRECISION    ZERO, S31, S32, S33, S34, DOTF
          DOTF(X,Y,S31)  =  DBLE(X)*DBLE(Y) + S31
          DATA      ZERO      /  0.D0  /
#endif
          DATA      HNAME               /  ' RFINV'  /
#include "finv.inc"
          RETURN
          END
