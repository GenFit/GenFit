*
* $Id: uginit.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: uginit.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:27  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.35  by  S.Giani
*-- Author :
      SUBROUTINE UGINIT
*
************************************************************************
*                                                                      *
*              To initialise GEANT3 program and read data cards        *
*                                                                      *
************************************************************************
*
#include "geant321/gckine.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcphys.inc"
*
*     -----------------------------------------------------------------
*
*             Initialize GEANT
C..geant..
      CALL G3INIT
*
*             Prints version number
*
      WRITE(LOUT,1000)
*
*             IKINE  = particle type (default=1=gamma)
*             PKINE(1)=particle energy
*             IKINE and PKINE can be changed with the data card KINE
*
      PKINE(1)=10.
      PKINE(5)=4.
      IKINE=1
      ICOMP=1
      IPAIR=1
      IBREM=1
      IANNI=1
      IPHOT=1
      IHADR=1
*
*             Initialize GEANT/ZBOOK data structures
*
C..geant..
      CALL G3ZINIT
*
*             Initialize graphics package
*
      CALL GDINIT
*
 1000 FORMAT(/,'  MODE VERSION 1.00 : ',/)
      END
