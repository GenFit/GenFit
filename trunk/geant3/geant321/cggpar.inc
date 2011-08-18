*
* $Id: cggpar.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cggpar.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGGPAR_INC
#define CERNLIB_GEANT321_CGGPAR_INC
*
*
* cggpar.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (LCGHEA=3, KCGSIZ=1)
      PARAMETER (          KCGATT=2)
      PARAMETER (          KCGNF =3)
      PARAMETER (LCGFAC=6, KCGAF =1)
      PARAMETER (          KCGAA =2)
      PARAMETER (          KCGBB =3)
      PARAMETER (          KCGCC =4)
      PARAMETER (          KCGDD =5)
      PARAMETER (          KCGNE =6)
      PARAMETER (LCGEDG=7, KCGAE =1)
      PARAMETER (          KCGX1 =2)
      PARAMETER (          KCGY1 =3)
      PARAMETER (          KCGZ1 =4)
      PARAMETER (          KCGX2 =5)
      PARAMETER (          KCGY2 =6)
      PARAMETER (          KCGZ2 =7)
#endif
*
* cggpar.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGGPAR.
#include "geant321/incend.inc"
*
#endif

#endif
