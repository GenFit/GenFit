*
* $Id: mzct.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mzct.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:53  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_MZCT_INC
#define CERNLIB_ZEBRA_MZCT_INC
*
*                 Memory Table parameters
*
* mzct.inc
*
      COMMON /MZCT/  MQDVGA,MQDVWI,JQSTMV,JQDVM1,JQDVM2,NQDVMV,IQFLIO
     +,              MQDVAC,NQNOOP,IQPART,NQFREE, IQTBIT,IQTVAL
     +,              IQTNMV,JQGAPM,JQGAPR,NQGAPN,NQGAP,IQGAP(5,4)
     +,              LQTA,LQTE, LQRTA,LQTC1,LQTC2,LQRTE
     +,                         LQMTA,LQMTB,LQMTE,LQMTLU,LQMTBR
     +,              LQMTC1,LQMTC2, NQFRTC,NQLIVE

#endif
