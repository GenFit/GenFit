*
* $Id: tzread2.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: tzread2.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:56  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_TZREAD2_INC
#define CERNLIB_ZEBRA_TZREAD2_INC
*
*
* tzread2.inc
*
#if !defined(CERNLIB_CRAZYERR)
      READ (LUNTQ,LINE(JFMTC1:JFMTC2),END=96) (LQ(J),J=LPUTA,LPUTE-1)
#else
      IF (NAME(5).NE.3)  THEN
          READ (LUNTQ,LINE(JFMTC1:JFMTC2),ERR=64,END=96)
     +                (LQ(J),J=LPUTA,LPUTE-1)
        ELSE
          READ (LUNTQ,LINE(JFMTC1:JFMTC2),ERR=64,END=96)
     +                (Q(J),J=LPUTA-8,LPUTE-9)
        ENDIF
   64 CONTINUE
#endif

#endif
