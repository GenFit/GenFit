*
* $Id: zvfwkc.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zvfwkc.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:53  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_ZVFWKC_INC
#define CERNLIB_ZEBRA_ZVFWKC_INC
*
*               ZVERIF working area
*
* zvfwkc.inc
*
                     CHARACTER CHWARN*10, TEXTHD*50
      COMMON /ZVFWKC/MASKST,IFLHD,NFATAL,NWARN, LFAILS,LFAILD,IXVSTO
     +,              JVMODE,IFLST,IFLLA,IFLDV,IFLBK
     +,              LSTOLO,LSTOHI, JDESLA,LOCAR,LOCARR,LOCARE,LATMP
     +,              LDIVLO,LDIVHI,MASKTO
     +,                     LBKLS,IDBK(2),NBKNN(4),LBKLN,LBKNX,IBKXX
     +,              NAMEPR(2), CHWARN,TEXTHD

#endif
