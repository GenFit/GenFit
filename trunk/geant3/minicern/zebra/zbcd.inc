*
* $Id: zbcd.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zbcd.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:52  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_ZBCD_INC
#define CERNLIB_ZEBRA_ZBCD_INC
*
*                 BCD char. in CETA order with A1 format
*
* zbcd.inc
*
      COMMON /ZBCD/  IQNUM2(11),IQLETT(26),IQNUM(10),   IQPLUS,IQMINS
     +,              IQSTAR,IQSLAS,IQOPEN,IQCLOS,IQDOLL,IQEQU, IQBLAN
     +,              IQCOMA,IQDOT, IQNUMB,IQAPO, IQEXCL,IQCOLO,IQQUOT
     +,              IQUNDE,IQCLSQ,IQAND, IQAT,  IQQUES,IQOPSQ,IQGREA
     +,              IQLESS,IQREVE,IQCIRC,IQSEMI,IQPERC,  IQLOWL(26)
     +,              IQCROP,IQVERT,IQCRCL,IQNOT, IQGRAV,  IQILEG
     +,              NQHOL0,NQHOLL(95)

#endif
