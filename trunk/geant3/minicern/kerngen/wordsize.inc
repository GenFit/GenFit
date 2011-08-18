*
* $Id: wordsize.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: wordsize.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:18  mclareni
* Kernlib
*
*
*
* wordsize.inc
*
*    # of bits/word
*    # of characters/word
#if defined(CERNLIB_QMIRTD)
      PARAMETER   (NBITPW=64, NCHAPW=8, NBITPC=8)
#else
      PARAMETER   (NBITPW=32)
      PARAMETER   (NCHAPW=4)
#endif
