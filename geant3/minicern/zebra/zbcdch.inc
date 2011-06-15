*
* $Id: zbcdch.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: zbcdch.inc,v $
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
#ifndef CERNLIB_ZEBRA_ZBCDCH_INC
#define CERNLIB_ZEBRA_ZBCDCH_INC
*
*               BCD char. in CETA order as CHARACTER
*
* zbcdch.inc
*
                      CHARACTER  CQALLC*96
      COMMON /ZBCDCH/ CQALLC
                      CHARACTER*1  CQLETT(96), CQNUM(10)
                      EQUIVALENCE (CQLETT(1),CQALLC(1:1))
                      EQUIVALENCE (CQNUM(1), CQALLC(27:27))

#endif
