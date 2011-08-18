*
* $Id: kerngt.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: kerngt.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:47  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:30  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/02/15 17:49:58  mclareni
* Kernlib
*
*
      SUBROUTINE KERNGT (LUN)
C
C CERN PROGLIB# Z001    KERNGT          .VERSION KERNLNX  1.00  930829
C ORIG.  16/03/89
C
C--   Store Pam title into /SLATE/ for printing by GSHOWF

      COMMON /SLATE/ MCH(40)


      CALL UCTOH1   (
     + 'KERNLNX  1.03  1995/07/18 13.54           '       
     +,  MCH, 40)
      CALL GSHOWF (LUN)
      RETURN
      END
