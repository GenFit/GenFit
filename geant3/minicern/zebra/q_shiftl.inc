*
* $Id: q_shiftl.inc 223 2007-12-10 15:12:03Z brun $
*
* $Log: q_shiftl.inc,v $
* Revision 1.1.1.1  2002/07/24 15:56:28  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1  1996/04/18 16:15:25  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_SHIFTL_INC
#define CERNLIB_ZEBRA_Q_SHIFTL_INC
*
*
* q_shiftl.inc
*
#if defined(CERNLIB_APOLLO)
      ISHFTL (IZW,NZB) = LSHFT (IZW, NZB)

#elif defined(CERNLIB_CDC)
      ISHFTL (IZW,NZB) = SHIFT (IZW,NZB) .AND. MASK(60-NZB)

#elif defined(CERNLIB_CRAY)
      ISHFTL (IZW,NZB) = SHIFTL (IZW, NZB)

#elif (defined(CERNLIB_MSDOS))&&(defined(CERNLIB_QF_F2C))
       ISHFTL(IZW,NZB) = LSHIFT(IZW,NZB)

#elif defined(CERNLIB_QF_SIEM)
      ISHFTL (IZW,NZB) = SHFTL (IZW, NZB)

#elif (defined(CERNLIB_LINUX))&&(!(defined(__G95__)||defined(CERNLIB_GFORTRAN)))
      ISHFTL (IZW,NZB) = LSHIFT (IZW,NZB)

#elif (defined(CERNLIB_SUN))&&(!defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
#elif (defined(CERNLIB_SUN))&&(defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)

#elif defined(CERNLIB_QMVDS)
C   ISHFTL / ISHFTR  are intrinsics of vdS

#elif 1
      ISHFTL (IZW,NZB) = ISHFT (IZW, NZB)

#endif

#endif
