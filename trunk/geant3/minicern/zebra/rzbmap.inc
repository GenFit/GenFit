*
* $Id: rzbmap.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: rzbmap.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:55  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_RZBMAP_INC
#define CERNLIB_ZEBRA_RZBMAP_INC
*
*
* rzbmap.inc
*
*
*     Dummy RZBMAP sequence. To be overridden from cradle
*
      PARAMETER(KDMAX=9)
      COMMON/RZBMAP/IDATA(KDMAX)

#endif
