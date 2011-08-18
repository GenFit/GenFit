*
* $Id: q_sbit.inc 223 2007-12-10 15:12:03Z brun $
*
* $Log: q_sbit.inc,v $
* Revision 1.1.1.1  2002/07/24 15:56:28  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.2  1996/04/18 16:15:23  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_Q_SBIT_INC
#define CERNLIB_ZEBRA_Q_SBIT_INC
*
*
* q_sbit.inc
*

#if defined(CERNLIB_APOLLO)
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHFT(1,IZP-1)) )
     +                             , LSHFT(AND(MZ,1),IZP-1) )

#elif defined(CERNLIB_CDC)
      MSBIT (MZ,IZW,IZP) =(IZW .AND. .NOT.SHIFT(MASK(1),IZP))
     +                  .OR. SHIFT(MZ .AND. .NOT.MASK(59),IZP-1)

#elif defined(CERNLIB_CRAY)
      MSBIT (MZ,IZW,IZP)     = (IZW .AND. SHIFT(MASK(65),IZP))
     +                     .OR. SHIFT(MZ.AND.MASK(127),IZP-1)

#elif (defined(CERNLIB_CONVEX))&&(!defined(CERNLIB_DOUBLE))

#elif (defined(CERNLIB_CONVEX))&&(defined(CERNLIB_DOUBLE))

#elif defined(CERNLIB_MSDOS)

#if defined(CERNLIB_QF_NDP)
      MSBIT (MZZ,IZZW,IZZP) = IOR (IAND (IZZW, NOT(ISHFT(1,IZZP-1)) )
     +                              , ISHFT(IAND(MZZ,1),IZZP-1))
#elif defined(CERNLIB_QF_F2C)
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHIFT(1,IZP-1)) )
     +                              , LSHIFT(AND(MZ,1),IZP-1))
#endif

#elif defined(CERNLIB_QF_SIEM)
      MSBIT (MZ,IZW,IZP) = LOR (LAND (IZW, LCOMPL(SHFTL(1,IZP-1)) )
     +                                 , SHFTL(LAND(MZ,1),IZP-1) )

#elif (defined(CERNLIB_LINUX))&&(!(defined(__G95__)||defined(CERNLIB_GFORTRAN)))
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHIFT(1,IZP-1)))
     +                              ,LSHIFT(AND(MZ,1),IZP-1))

#elif defined(CERNLIB_PDP)
      MSBIT (MZ,IZW,IZP) = IOR (IAND (IZW, NOT(ISHFT(1,IZP-1)) )
     +                             , ISHFT (IAND(MZ,1),IZP-1) )

#elif defined(CERNLIB_SUN)
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHIFT(1,IZP-1)) )
     +                             , LSHIFT(AND(MZ,1),IZP-1) )

#elif defined(CERNLIB_VAXVMS)

#elif defined(CERNLIB_QMVDS)
      MSBIT (MZ,IZW,IZP) = OR (AND (IZW, NOT(LSHFT(1,IZP-1)) )
     +                             , LSHFT(AND(MZ,1),IZP-1) )

#elif 1
      MSBIT (MZ,IZW,IZP) = IOR (IAND (IZW, NOT(ISHFT(1,IZP-1)) )
     +                              , ISHFT(IAND(MZ,1),IZP-1))
#endif

#endif
