*
* $Id: cgxcom.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cgxcom.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_CGXCOM_INC
#define CERNLIB_GEANT321_CGXCOM_INC
*
*
* cgxcom.inc
*
#if !defined(CERNLIB_INCLUDE)
      COMMON /CGXCOM/ IXLUN,NXMAX,NXCUR,NXNODE,NXEDGE,NXFACE
#endif
*
* cgxcom.inc
*
#if defined(CERNLIB_INCLUDE)
#include "geant321/incbeg.inc"
       +   CGXCOM.
#include "geant321/incend.inc"
*
#endif

#endif
