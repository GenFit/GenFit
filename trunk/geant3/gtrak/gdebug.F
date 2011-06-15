*
* $Id: gdebug.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gdebug.F,v $
* Revision 1.2  2003/11/28 11:23:56  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:41  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:40  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE G3DEBUG
*
************************************************************************
*             Standard GEANT debug                                     *
*                                                                      *
*          ==>Called by : <USER>,GUSTEP                                *
*             Author    R.Brun  *********                              *
************************************************************************
*
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
*
*             Debug/plot event
      IF(IDEBUG.NE.0) THEN
        IF((ISWIT(2).EQ.1).OR.(ISWIT(3).EQ.1)
     +                    .OR.(ISWIT(2).EQ.4)) CALL G3SXYZ
        IF (ISWIT(2).EQ.2) CALL G3PCXYZ
        IF (ISWIT(1).EQ.2) CALL G3PGKIN
        IF (ISWIT(2).EQ.3) THEN
           IF(ISWIT(4).EQ.3.AND.CHARGE.EQ.0.)RETURN
           CALL GDCXYZ
*           CALL IGSA(0)
        ELSEIF(ISWIT(2).EQ.4) THEN
           IF(ISTOP.NE.0) CALL GDTRAK(' ')
        ENDIF
      ENDIF
      END
 
