*
* $Id: q_andor.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: q_andor.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.2  1997/02/04 17:35:38  mclareni
* Merge Winnt and 97a versions
*
* Revision 1.1.1.1.2.1  1997/01/21 11:30:28  mclareni
* All mods for Winnt 96a on winnt branch
*
* Revision 1.1.1.1  1996/02/15 17:49:19  mclareni
* Kernlib
*
*
*    IAND / IOR / IEOR / NOT  standard?
*
* q_andor.inc
*
#if defined(CERNLIB_QMSUN)||defined(CERNLIB_QMLNX)
      IAND (IZV,IZW) = AND (IZV, IZW)
      IOR  (IZV,IZW) =  OR (IZV, IZW)
      IEOR (IZV,IZW) = XOR (IZV, IZW)
#elif defined(CERNLIB_QMDOS)&&defined(CERNLIB_QF2C)
      IAND(IZV,IZW)  = AND(IZV,IZW)
       IOR(IZV,IZW)  =  OR(IZV,IZW)
      IEOR(IZV,IZW)  = XOR(IZV,IZW)
#else
#endif
