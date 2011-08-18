*
* $Id: zceta.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zceta.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_ZCETA_INC
#define CERNLIB_ZEBRA_ZCETA_INC
*
*                CETA values in internal order
*
* zceta.inc
*
#if defined (CERNLIB_CDC)
      PARAMETER     (NQTCET=64)
      COMMON /ZCETA/ IQCETA(64), IQTCET(64)

#elif defined (CERNLIB_PDP)
      PARAMETER     (NQTCET=128)
      COMMON /ZCETA/ IQCETA(128),IQTCET(128)

#elif defined (CERNLIB_QMUNO)
      PARAMETER     (NQTCET=512)
      COMMON /ZCETA/ IQCETA(512),IQTCET(512)

#elif 1
      PARAMETER     (NQTCET=256)
      COMMON /ZCETA/ IQCETA(256),IQTCET(256)

#endif

#endif
