*
* $Id: q_jbit.inc 221 2007-11-25 08:20:46Z brun $
*
* $Log: q_jbit.inc,v $
* Revision 1.2  2002/12/02 16:37:45  brun
* Changes from Federico Carminati and Peter Hristov who ported the system
* on the Ithanium processors.It is tested on HP, Sun, and Alpha, everything
* seems to work. The optimisation is switched off in case of gcc2.xx.yyy
*
* Revision 1.1.1.1  2002/07/24 15:56:28  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.2  1996/04/18 16:15:06  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_JBIT_INC
#define CERNLIB_ZEBRA_Q_JBIT_INC
*
*
* q_jbit.inc
*
#if defined(CERNLIB_ALLIANT)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_APOLLO)
      JBIT (IZW,IZP)    =   AND (RSHFT(IZW,IZP-1), 1)

#elif defined(CERNLIB_CDC)
      JBIT (IZW,IZP)     = SHIFT(IZW,61-IZP) .AND. .NOT.MASK(59)

#elif defined(CERNLIB_CRAY)
      JBIT (IZW,IZP)     = SHIFTR(IZW,IZP-1) .AND. MASK(127)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      JBIT (IZZW,IZZP)      = IBITS (IZZW,IZZP-1,1)
#elif defined(CERNLIB_QF_F2C)
      JBIT (IZW,IZP)      = AND(RSHIFT(IZW,IZP-1),1)
#endif

#elif defined(CERNLIB_HPUX)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_QF_SIEM)
      JBIT (IZW,IZP)    = SHFTR (SHFTL(IZW,32-IZP), 31)

#elif (defined(CERNLIB_LINUX))&&(!defined(__G95__))&&(!defined(CERNLIB_GFORTRAN))
      JBIT(IZW,IZP)     = AND(RSHIFT(IZW,IZP-1),1)

#elif defined(CERNLIB_OS9)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif defined(CERNLIB_PDP)
      JBIT (IZW,IZP)     = IBITS (IZW,IZP-1,1)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      JBIT (IZW,IZP)    = AND (lrshft(IZW,IZP-1), 1)

#elif defined(CERNLIB_QMUNO)
      JBIT (IZW,IZP)     = BITS (IZW,37-IZP,1)

#elif defined(CERNLIB_QMVDS)
      JBIT (IZW,IZP)  =  JBYT (IZW,IZP,1)

#elif defined(CERNLIB_VAXVMS)

#elif 1
      JBIT (IZW,IZP)     =  IAND (ISHFT(IZW,1-IZP), 1)

#endif

#endif
