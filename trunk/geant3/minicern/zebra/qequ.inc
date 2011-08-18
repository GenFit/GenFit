*
* $Id: qequ.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: qequ.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:27  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:53  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_QEQU_INC
#define CERNLIB_ZEBRA_QEQU_INC
*
*
* qequ.inc
*
              DIMENSION  LQLUP(99), LQLORG(98),IQWIDN(97),IQWIDH(96)
     +,                  IQWNL(95), IQWNS(94), IQWND(93)
            EQUIVALENCE (LQLUP(8),  LQLORG(7), IQWIDN(6), IQWIDH(5)
     +,                  IQWNL(4),  IQWNS(3),  IQWND(2), IQ(1))

#endif
