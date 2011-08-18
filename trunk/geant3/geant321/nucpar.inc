*
* $Id: nucpar.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: nucpar.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:41  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_NUCPAR_INC
#define CERNLIB_GEANT321_NUCPAR_INC
*
*
* nucpar.inc
*
*
*=== nucpar ===========================================================*
*
*
*----------------------------------------------------------------------*
*     Created on 6-9-89 by A. Ferrari                                  *
*                       included in:                                   *
*                                    nucevt                            *
*                                    evxtes                            *
*                                    distnu                            *
*                                    eventq                            *
*----------------------------------------------------------------------*
 
      PARAMETER (MXPNUC=600)
      CHARACTER*8 ANNU
      COMMON /FKNUPA/  PXNU(MXPNUC),   PYNU(MXPNUC),   PZNU(MXPNUC),
     &                HEPNU(MXPNUC),   AMNU(MXPNUC),   ANNU(MXPNUC),
     &                ICHNU(MXPNUC), IBARNU(MXPNUC),  NRENU(MXPNUC)
 
 

#endif
