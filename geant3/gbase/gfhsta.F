*
* $Id: gfhsta.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gfhsta.F,v $
* Revision 1.2  2003/11/28 11:23:55  brun
* New version of geant321 with all geant3 routines renamed from G to G3
*
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.2  1997/10/24 16:17:42  mclareni
* NZLEFT 2nd parameter changed to 0 from blank
*
* Revision 1.1.1.1  1995/10/24 10:20:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3FHSTA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Fills histograms statistics                              *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUOUT                                *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#ifdef NEVER
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gclist.inc"
      COMMON/QUEST/IQUEST(100)
C.
C.    ------------------------------------------------------------------
C.
      IF(NHSTA.GT.0) THEN
         CALL GLOOK('TIME',LHSTA,NHSTA,ID)
         IF(ID.NE.0)THEN
            CALL TIMED(TEVENT)
            CALL HFILL(-ID,TEVENT,0.,1.)
         ENDIF
*
         CALL GLOOK('SIZE',LHSTA,NHSTA,ID)
         IF(ID.NE.0)THEN
            NL=NZLEFT(IXDIV,0)
            CALL HFILL(-ID,IQUEST(11)+0.5,0.,1.)
         ENDIF
*
         CALL GLOOK('MULT',LHSTA,NHSTA,ID)
         IF(ID.NE.0)THEN
            CALL HFILL(-ID,NTMULT+0.5,0.,1.)
         ENDIF
*
         CALL GLOOK('NTRA',LHSTA,NHSTA,ID)
         IF(ID.NE.0)THEN
            CALL HFILL(-ID,NTRACK+0.5,0.,1.)
         ENDIF
*
         CALL GLOOK('STAK',LHSTA,NHSTA,ID)
         IF(ID.NE.0)THEN
            CALL HFILL(-ID,NSTMAX+0.5,0.,1.)
         ENDIF
*
      ENDIF
#endif
      END
