*
* $Id: q_sbit1.inc 223 2007-12-10 15:12:03Z brun $
*
* $Log: q_sbit1.inc,v $
* Revision 1.1.1.1  2002/07/24 15:56:28  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1  1996/04/18 16:15:24  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_SBIT1_INC
#define CERNLIB_ZEBRA_Q_SBIT1_INC
*
*
* q_sbit1.inc
*

#if defined(CERNLIB_APOLLO)
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHFT(1,IZP-1)  )

#elif defined(CERNLIB_CDC)
      MSBIT1 (IZW,IZP)   = IZW .OR.       SHIFT(MASK(1),IZP)

#elif defined(CERNLIB_CRAY)
      MSBIT1 (IZW,IZP)       =  IZW .OR.  SHIFT(MASK(1), IZP)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))
      MSBIT1 (IZW,IZP)   = IBSET (IZW,IZP-1)

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))
      MSBIT1 (IZW,IZP)   = IBSET (IZW,IZP-1)

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      MSBIT1 (IZZW,IZZP)   = IOR  (IZZW,     ISHFT(1,IZZP-1)  )
#elif defined(CERNLIB_QF_F2C)
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHIFT(1,IZP-1)  )
#endif

#elif defined(CERNLIB_QF_SIEM)
      MSBIT1 (IZW,IZP)   = LOR  (IZW,        SHFTL(1,IZP-1)  )

#elif (defined(CERNLIB_LINUX))&&(!(defined(__G95__)||defined(CERNLIB_GFORTRAN)))
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHIFT(1,IZP-1))

#elif defined(CERNLIB_PDP)
      MSBIT1 (IZW,IZP)   = IOR  (IZW,     ISHFT(1,IZP-1) )

#elif defined(CERNLIB_SUN)
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHIFT(1,IZP-1)  )

#elif defined(CERNLIB_VAXVMS)

#elif defined(CERNLIB_QMVDS)
      MSBIT1 (IZW,IZP)   = OR  (IZW,     LSHFT(1,IZP-1)  )

#elif 1
      MSBIT1 (IZW,IZP)   = IOR  (IZW,     ISHFT(1,IZP-1)  )
#endif

#endif
