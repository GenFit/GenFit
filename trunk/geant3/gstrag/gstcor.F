*
* $Id: gstcor.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gstcor.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:39  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE G3STCOR(CC,POT,X,DENS,DELTA)
 
C.
C.    *****************************************************************
C.    *                                                               *
C.    *      Sternheimer correction as computed in G3DRELX            *
C.    *                                                               *
C.    *    ==>Called by :                                             *
C.    *       Author    L.Urban    *********                          *
C.    *                                                               *
C.    *****************************************************************
*
*         condensed material ? ( dens .gt. 0.05 ? )
      IF(DENS.GT.0.05) THEN
*        IF(POT.LT.1.E-7) THEN   (POT was in GeV, now in keV)
         IF(POT.LT.0.1) THEN
            IF(CC.LT.3.681) THEN
               X0=0.2
            ELSE
               X0=0.326*CC-1.
            ENDIF
            X1=2.
         ELSE
            IF(CC.LT.5.215) THEN
               X0=0.2
            ELSE
               X0=0.326*CC-1.5
            ENDIF
            X1=3.
         ENDIF
*         gas ?   ( dens . le . 0.05 ? )
      ELSE
         IF(CC.LE.12.25) THEN
            IP=INT((CC-10.)/0.5)+1
            IF(IP.LT.0) IP=0
            IF(IP.GT.4) IP=4
            X0=1.6+0.1*FLOAT(IP)
            X1=4.
         ELSE
            IF(CC.LE.13.804) THEN
               X0=2.
               X1=5.
            ELSE
               X0=0.326*CC-2.5
               X1=5.
            ENDIF
         ENDIF
      ENDIF
*
      XA=CC/4.606
      XM=3.
      AA=4.606*(XA-X0)/(X1-X0)**XM
*
      DELTA=0.
      IF(X.GT.X0) THEN
         DELTA=4.606*X-CC
         IF(X.LT.X1) DELTA=DELTA+AA*(X1-X)**XM
      ENDIF
 
      END
