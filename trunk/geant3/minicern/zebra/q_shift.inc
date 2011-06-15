*
* $Id: q_shift.inc 223 2007-12-10 15:12:03Z brun $
*
* $Log: q_shift.inc,v $
* Revision 1.2  2002/12/02 16:37:45  brun
* Changes from Federico Carminati and Peter Hristov who ported the system
* on the Ithanium processors.It is tested on HP, Sun, and Alpha, everything
* seems to work. The optimisation is switched off in case of gcc2.xx.yyy
*
* Revision 1.1.1.1  2002/07/24 15:56:28  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_SHIFT_INC
#define CERNLIB_ZEBRA_Q_SHIFT_INC
*
*
* q_shift.inc
*
#if defined(CERNLIB_APOLLO)
      ISHFTL (IZW,NZB) = LSHFT (IZW, NZB)
      ISHFTR (IZW,NZB) = RSHFT (IZW, NZB)

#elif defined(CERNLIB_CDC)
      ISHFTL (IZW,NZB) = SHIFT (IZW,NZB) .AND. MASK(60-NZB)
      ISHFTR (IZW,NZB) = SHIFT (IZW.AND.MASK(60-NZB), 60-NZB)

#elif defined(CERNLIB_CRAY)
      ISHFTL (IZW,NZB) = SHIFTL (IZW, NZB)
      ISHFTR (IZW,NZB) = SHIFTR (IZW, NZB)

#elif (defined(CERNLIB_MSDOS))&&(defined(CERNLIB_QF_F2C))
C--
C--    Note: RSHIFT is arithmetic with F2C, it cannot be used
C
       ISHFTL(IZW,NZB) = LSHIFT(IZW,NZB)
C       ISHFTR(IZW,NZB) = RSHIFT(IZW,NZB)
C       ISHFT(IZW,NZB)  = LSHIFT(IZW,NZB)

#elif defined(CERNLIB_QF_SIEM)
      ISHFTL (IZW,NZB) = SHFTL (IZW, NZB)
      ISHFTR (IZW,NZB) = SHFTR (IZW, NZB)

#elif (defined(CERNLIB_LINUX))&&(!(defined(__G95__)||defined(CERNLIB_GFORTRAN)))
      ISHFTL (IZW,NZB) = LSHIFT (IZW,NZB)
      ISHFTR (IZW,NZB) = RSHIFT (IZW,NZB)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
      ISHFTR (IZW,NZB) = lrshft (IZW, NZB)
#elif (defined(CERNLIB_SUN))&&(defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
      ISHFTR (IZW,NZB) =  ishft (IZW,-NZB)

#elif defined(CERNLIB_QMVDS)
C   ISHFTL / ISHFTR  are intrinsics of vdS

#elif 1
      ISHFTL (IZW,NZB) = ISHFT (IZW, NZB)
      ISHFTR (IZW,NZB) = ISHFT (IZW,-NZB)

#endif

#endif
