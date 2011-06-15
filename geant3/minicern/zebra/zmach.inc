*
* $Id: zmach.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zmach.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_ZMACH_INC
#define CERNLIB_ZEBRA_ZMACH_INC
*
*                Machine parameters
*
* zmach.inc
*
#include "zebra/zmachfix.inc"
      COMMON /ZMACH/ NQBITW,NQBITC,NQCHAW
     +,              NQLNOR,NQLMAX,NQLPTH,NQRMAX,IQLPCT,IQNIL

#endif
