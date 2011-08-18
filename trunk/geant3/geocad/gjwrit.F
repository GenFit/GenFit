*
* $Id: gjwrit.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gjwrit.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GJWRIT
C
C    Writes 80 characters in a line and skips spaces
C
C    CALLED BY ALL SHAPE ROUTINES in CADINT
C
C    BLKSTR = The character string to be written
C    NPLACE = The place of the current character in the line
C
C---------------------------------------------------------------------
C
#include "geant321/gcsetf.inc"
 
C
C---------------------------------------------------------------------
C
C
      J=NPLACE
C
      LBLK=LEN(BLKSTR)
      DO 10  I=1, LBLK
         IF (BLKSTR(I:I) .NE. ' ') THEN
            WRITE(JVSTR(J:J),'(A1)')BLKSTR(I:I)
            IF (J .EQ. 80) THEN
               WRITE(NUNIT1,'(A80)')JVSTR
               JVSTR=' '
               J=0
            END IF
            J=J+1
         END IF
   10 CONTINUE
C
      NPLACE=J
C
      BLKSTR=' '
C
      END
