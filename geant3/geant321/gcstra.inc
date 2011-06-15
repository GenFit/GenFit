*
* $Id: gcstra.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gcstra.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:33  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_GCSTRA_INC
#define CERNLIB_GEANT321_GCSTRA_INC
*
*
* gcstra.inc
*
#if !defined(CERNLIB_INCLUDE)
      PARAMETER (ILTAB=200)
      COMMON /GCSTRA / EMAX,EM(200),SFINT,EPSR(ILTAB),EPSI(ILTAB),
     +  FINT(ILTAB),EMIN,EPPS,BETA2,GAMMA2,WP2,S2,MEEV,EMM(200),
     +  GAMLOG(21),NP,NTAB,IE,NFACT,NICOLL
*
      INTEGER NP,NTAB,IE,NFACT
      REAL    EMAX,EM,SFINT,EPSR,EPSI,FINT,EMIN,EPPS,EMM,GAMLOG
     +        ,BETA2,GAMMA2,WP2,S2,MEEV
*
#endif

#endif
