*
* $Id: aadat.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: aadat.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_AADAT_INC
#define CERNLIB_GEANT321_AADAT_INC
*
*
* aadat.inc
*
*                                                                      *
*=== aadat ============================================================*
*                                                                      *
      PARAMETER (IPROM = 100)
      PARAMETER (ITARM = 100)
      COMMON/FKAADA/ENPERN,PLPERN,SIGAA(IPROM,ITARM),SIGNN,TMASS,AAEVNO,
     +             SELAA(IPROM,ITARM),RLASTP(IPROM,ITARM),
     +             MATPRO,MATTAR,LASTM,IPROA,IPROZ,ITARA,ITARZ
 

#endif
