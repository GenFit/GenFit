*
* $Id: gstdn.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gstdn.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      FUNCTION G3STDN(LGE)
*   ********************************************************************
*   *                                                                  *
*   *    dN^2/dxdE function for the preparation of the table as in     *
*   *    Grishin,Ermilova,Kotelnikov,NIM A307(1991),273                *
*   *                                                                  *
*   *    ==> Called by : G3STTAB                                       *
*   *        Author    : P. Nevski                                     *
*   *                                                                  *
*   ********************************************************************
#include "geant321/gconsp.inc"
#include "geant321/gcstra.inc"
*
      REAL LGE
      COMPLEX I,CEPS1,C1,C2
      SEXP(X) = EXP(MAX(X,-99.))
*
*     ------------------------------------------------------------------
      I=(0.,1.)
      E=EXP(LGE)
*
*    Epsilon-1
      E0     = MAX(LGE,EMIN)
      E2INV  = 1/E**2
      CEPS1  =  CMPLX(G3STINT(LGE,EMM,EPSR,NTAB)*E2INV,
     +       SEXP(G3STINT(E0,EMM,EPSI,NTAB)))
      C1     = 1/GAMMA2-CEPS1*BETA2
      C2     = C1/(1+CEPS1)*LOG(2*BETA2*MEEV/(E*C1))
      G3STDN = S2/BETA2 * (-2*AIMAG(C2)/(WP2*PI) + (1-SEXP(G3STINT(E0
     +, EMM,FINT,NTAB)))*E2INV) * E
      END
