*
* $Id: bankparq.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: bankparq.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:49  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:56  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_BANKPARQ_INC
#define CERNLIB_ZEBRA_BANKPARQ_INC
*
*  BANK  PARAMETERS
*
* bankparq.inc
*
      PARAMETER(NBKHDQ=1,NBKCTQ=9,NBKOHQ=NBKHDQ+NBKCTQ,NBKJMQ=12)
      PARAMETER(MBKLNQ=0,JBKSKQ=1,NBKSKQ=16)
      PARAMETER(MBKIOQ=MBKLNQ,JBKIOQ=NBKSKQ+1,NBKIOQ=16)
      PARAMETER(JBKUSQ=1,NBKUSQ=18,JBKSYQ=NBKUSQ+1,NBKSYQ=32-NBKUSQ)

#endif
