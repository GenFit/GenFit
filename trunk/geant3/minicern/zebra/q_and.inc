*
* $Id: q_and.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: q_and.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1  1996/04/18 16:15:05  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_AND_INC
#define CERNLIB_ZEBRA_Q_AND_INC
*
* q_and.inc
*
#if defined(CERNLIB_APOLLO)
      IAND (IZV,IZW) = AND (IZV, IZW)

#elif defined(CERNLIB_CDC)
      IAND (IZV,IZW) = IZV .AND. IZW

#elif defined(CERNLIB_MSDOS)&&defined(CERNLIB_QF_F2C)
      IAND(IZV,IZW)  = AND(IZV,IZW)

#elif defined(CERNLIB_QMIBMSI)
      IAND (IZV,IZW) = LAND (IZV, IZW)

#elif defined(CERNLIB_LNX)
      IAND(IZV,IZW)  = AND(IZV,IZW)

#elif defined(CERNLIB_SUN)
      IAND (IZV,IZW) = AND (IZV, IZW)
 
#elif 1

#endif

#endif
