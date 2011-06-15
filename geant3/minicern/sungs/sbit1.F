*
* $Id: sbit1.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: sbit1.F,v $
* Revision 1.1.1.1  1999/05/18 15:55:32  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      SUBROUTINE SBIT1 (IZW,IZP)
C
C CERN PROGLIB# M421    SBIT1           .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
      IZW = OR (IZW, LSHIFT(1,IZP-1))
      RETURN
      END
