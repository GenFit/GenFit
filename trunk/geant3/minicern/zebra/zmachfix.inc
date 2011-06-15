*
* $Id: zmachfix.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zmachfix.inc,v $
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
#ifndef CERNLIB_ZEBRA_ZMACHFIX_INC
#define CERNLIB_ZEBRA_ZMACHFIX_INC
*
*
* zmachfix.inc
*
#if defined (CERNLIB_UNIX)
      PARAMETER      (IQBITW=32, IQBITC=8, IQCHAW=4)

#elif defined (CERNLIB_CDC)
      PARAMETER      (IQBITW=60, IQBITC=6, IQCHAW=10)

#elif defined (CERNLIB_CRAY)
      PARAMETER      (IQBITW=64, IQBITC=8, IQCHAW=8)

#elif (defined (CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      PARAMETER      (IQBITW=64, IQBITC=8, IQCHAW=8)

#elif defined (CERNLIB_PDP)
      PARAMETER      (IQBITW=36, IQBITC=7, IQCHAW=5)

#elif defined (CERNLIB_QMUNO)
      PARAMETER      (IQBITW=36, IQBITC=9, IQCHAW=4)

#elif 1
      PARAMETER      (IQBITW=32, IQBITC=8, IQCHAW=4)

#endif

#endif
