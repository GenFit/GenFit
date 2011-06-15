*
* $Id: gxclos.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gxclos.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXCLOS(LUN)
*
*                To close unit LUN
*
#if defined(CERNLIB_IBM)
      CHARACTER*80 CHTEMP
      CHARACTER*3 CHLUN
#endif
*________________________________________________________
*
         CLOSE (LUN)
#if defined(CERNLIB_IBM)
         WRITE(CHLUN,'(I3)')LUN
         CHTEMP='FILEDEF '//CHLUN//'  CLEAR'
         CALL VMCMS(CHTEMP,ISTAT)
#endif
*
      END
 
