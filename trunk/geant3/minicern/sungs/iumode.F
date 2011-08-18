*
* $Id: iumode.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: iumode.F,v $
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      FUNCTION IUMODE (IWORD)
C
C CERN PROGLIB# M506    IUMODE          .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
C
      J = RSHIFT (IWORD,23)
      IF (J.EQ.-1)  J=0
C
      IUMODE= J
      RETURN
      END
