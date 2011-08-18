*
* $Id: intrac.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: intrac.F,v $
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:39  mclareni
* Kernlib
*
*
      LOGICAL FUNCTION INTRAC(DUMMY)
C
C CERN PROGLIB# Z044    INTRAC          .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
C-    The program is considered interactive if standard input
C-    is connected to the terminal

      LOGICAL      ISATTY

      INTRAC = ISATTY (5)
      RETURN
      END
