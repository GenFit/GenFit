*
* $Id: chkzeb.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: chkzeb.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:43  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE CHKZEB(NW,IX)
C
C Check if NW words are available in ZEBRA division IX
C
C ZEBRA user communication common
      COMMON/ QUEST / IQUEST(100)
C
      CALL MZNEED(IX,NW,'G')
      IF(IQUEST(11).LT.0) THEN
         PRINT *,'******************************************'
         PRINT *,'*            G C A L O R                 *'
         PRINT *,'*   NOT enough space available in ZEBRA  *'
         PRINT '('' *  division '',I3,'' to store '',I8,               '
     +   //'            '' words  *'')',IX,NW
         PRINT *,'*                                        *'
         PRINT *,'*  INCREASE ZEBRA COMMON SIZE AND RERUN  *'
         PRINT *,'*                                        *'
         PRINT *,'*             RUN TERMINATED             *'
         PRINT *,'******************************************'
         STOP
      ENDIF
      RETURN
      END
