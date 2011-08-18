*
* $Id: gfcoul.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfcoul.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      FUNCTION G3FCOUL(Z)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Compute Coulomb correction for pair production and Brem       *
C.    *  REFERENCE : EGS MANUAL SLAC 210 - UC32 - JUNE 78              *
C.    *                        FORMULA 2.7.17                          *
C.    *                                                                *
C.    *    ==>Called by : GSMIXT                                       *
C.    *       Author    M.Maire *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
      PARAMETER (ALPHA = 7.29927E-03)
C
      AZ2   = (ALPHA*Z)**2
      AZ4   =   AZ2 * AZ2
      FP    = ( 0.0083*AZ4 + 0.20206 + 1./(1.+AZ2) ) * AZ2
      FM    = ( 0.0020*AZ4 + 0.0369  ) * AZ4
      G3FCOUL = FP - FM
      END
