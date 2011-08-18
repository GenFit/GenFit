*
* $Id: q_xor.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: q_xor.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1  1996/04/18 16:15:26  mclareni
* Incorporate changes from J.Zoll for version 3.77
*
*
*
#ifndef CERNLIB_ZEBRA_Q_XOR_INC
#define CERNLIB_ZEBRA_Q_XOR_INC
*
* q_xor.inc
*
#if defined(CERNLIB_APOLLO)
      IEOR (IZV,IZW) = XOR (IZV, IZW)

#elif defined(CERNLIB_CDC)
      IEOR (IZV,IZW) = IZV .XOR. IZW

#elif defined(CERNLIB_MSDOS)&&defined(CERNLIB_QF_F2C)
      IEOR(IZV,IZW)  = XOR(IZV,IZW)

#elif defined(CERNLIB_QMIBMSI)
      IEOR (IZV,IZW) = LXOR (IZV, IZW)

#elif defined(CERNLIB_LNX)
      IEOR(IZV,IZW)  = XOR(IZV,IZW)

#elif defined(CERNLIB_SUN)
      IEOR (IZV,IZW) = XOR (IZV, IZW)
 
#elif 1

#endif

#endif
