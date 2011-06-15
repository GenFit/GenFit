*
* $Id: mzca.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mzca.inc,v $
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
#ifndef CERNLIB_ZEBRA_MZCA_INC
#define CERNLIB_ZEBRA_MZCA_INC
*
*     Process Master parameters
*
* mzca.inc
*
#if defined(CERNLIB_QTRHOLL)
      COMMON /MZCA/  NQSTOR,NQOFFT(16),NQOFFS(16),NQALLO(16), NQIAM
     +,              LQATAB,LQASTO,LQBTIS, LQWKTB,NQWKTB,LQWKFZ
     +,              MQKEYS(3),NQINIT,NQTSYS,NQM99,NQPERM,NQFATA,NQCASE
     +,              NQTRAC,MQTRAC(48)
                                       EQUIVALENCE (KQSP,NQOFFS(1))
#endif
*    Process Master parameters
#if !defined(CERNLIB_QTRHOLL)
      COMMON /MZCA/  NQSTOR,NQOFFT(16),NQOFFS(16),NQALLO(16), NQIAM
     +,              LQATAB,LQASTO,LQBTIS, LQWKTB,NQWKTB,LQWKFZ
     +,              MQKEYS(3),NQINIT,NQTSYS,NQM99,NQPERM
     +,              NQFATA,NQCASE,NQTRAC
                                       EQUIVALENCE (KQSP,NQOFFS(1))
      COMMON /MZCA2/ MQTRAC(44)
                     CHARACTER  MQTRAC*8
#endif

#endif
