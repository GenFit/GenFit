*
* $Id: q_shift.inc 223 2007-12-10 15:12:03Z brun $
*
* $Log: q_shift.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.2  1996/03/05 17:04:19  cernlib
* Replace #else by elif 1 ( for SGI cpp )
*
* Revision 1.1.1.1  1996/02/15 17:49:19  mclareni
* Kernlib
*
*
*
* q_shift.inc
*
#if defined(CERNLIB_QMSUN)&&defined(CERNLIB_BUGLRSHFT)
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
      ISHFTR (IZW,NZB) = ishft (IZW, -NZB)
#elif defined(CERNLIB_QMSUN)&&(!defined(CERNLIB_BUGLRSHFT))
      ISHFTL (IZW,NZB) = LSHIFT (IZW, NZB)
      ISHFTR (IZW,NZB) = lrshft (IZW, NZB)
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF_F2C)
C--
C--    Note: F2C-ish RSHIFT is arithmetic right shift.
C--          This means sign bits are shifted in at the high end
C
       ISHFTL(IZW,NZB) = LSHIFT(IZW,NZB)
C       ISHFTR(IZW,NZB) = RSHIFT(IZW,NZB)
C       ISHFT(IZW,NZB)  = LSHIFT(IZW,NZB)
#elif (defined(CERNLIB_QMLNX)&&!defined(CERNLIB_GFORTRAN))
      ISHFTL (IZW,NZB) = LSHIFT (IZW,NZB)
#elif 1
      ISHFTL (IZW,NZB) = ISHFT (IZW,  NZB)
      ISHFTR (IZW,NZB) = ISHFT (IZW, -NZB)
#endif
