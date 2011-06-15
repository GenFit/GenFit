*
* $Id: jbit.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: jbit.F,v $
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      FUNCTION JBIT (IZW,IZP)
C
C CERN PROGLIB# M421    JBIT            .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
      JBIT = AND (RSHIFT(IZW,IZP-1), 1)
      RETURN
      END
