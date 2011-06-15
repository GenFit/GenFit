*
* $Id: rzcstr.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: rzcstr.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:55  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_RZCSTR_INC
#define CERNLIB_ZEBRA_RZCSTR_INC
*
*
* rzcstr.inc
*
      PARAMETER     (MAXFILES=128, MAXSTRIP=21)
      CHARACTER*128  RZNAMES(MAXFILES),RZSFILE(MAXSTRIP)
      COMMON/RZCSTRC/RZNAMES,RZSFILE
      COMMON/RZCSTRI/ISLAST,ISTRIP(MAXFILES),NSTRIP(MAXFILES),
     +                      NRSTRIP(MAXFILES)
C

#endif
