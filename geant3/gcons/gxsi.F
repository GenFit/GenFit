*
* $Id: gxsi.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gxsi.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      FUNCTION G3XSI (Z)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  COMPUTE SCREENING FACTOR   FOR PAIR PRODUCTION AND BREM       *
C.    *  REFERENCE : EGS MANUAL SLAC 210 - UC32 - JUNE 78              *
C.    *                        FORMULA 2.7.22                          *
C.    *                                                                *
C.    *    ==>Called by : GSMIXT                                       *
C.    *       Author    M.Maire *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      DATA AL183  , AL1440 / 5.20948 , 7.27239 /
      ALZ  = LOG(Z)/3.
      G3XSI  = (AL1440 - 2*ALZ) / (AL183 - ALZ - G3FCOUL(Z))
      END
