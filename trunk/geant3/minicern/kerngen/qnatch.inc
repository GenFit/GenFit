*
* $Id: qnatch.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: qnatch.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:29  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:14  mclareni
* Kernlib
*
*
*  Siemens compiler
*
* qnatch.inc
*
#if defined(CERNLIB_QMIBM)
#if defined(CERNLIB_QF_SIEM)
      COMMON /QNATCH/ NATCH(256)
      INTEGER*2       NATCH
#else
      INTEGER*2       NATCH(256)
#endif
#include "qnatchdata.inc"
#endif
