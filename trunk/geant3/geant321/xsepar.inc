*
* $Id: xsepar.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: xsepar.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_XSEPAR_INC
#define CERNLIB_GEANT321_XSEPAR_INC
*
*
* xsepar.inc
*
*
*=== xsepar ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Xsec Parameters for neutrons and protons:                        *
*                                                                      *
*     Created on  20 september 1991  by  Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 20-sep-91     by    Alfredo Ferrari               *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                        BERRTP                                        *
*                        XSINEU                                        *
*                        XSIPRO                                        *
*                        PREPRE                                        *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON / FKXSEP / AANXSE (100), BBNXSE (100), CCNXSE (100),
     &                  DDNXSE (100), EENXSE (100), ZZNXSE (100),
     &                  EMNXSE (100), XMNXSE (100),
     &                  AAPXSE (100), BBPXSE (100), CCPXSE (100),
     &                  DDPXSE (100), EEPXSE (100), FFPXSE (100),
     &                  ZZPXSE (100), EMPXSE (100), XMPXSE (100)
 

#endif
