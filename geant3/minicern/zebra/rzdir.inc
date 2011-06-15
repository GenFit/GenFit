*
* $Id: rzdir.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: rzdir.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:50  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:28  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 10:46:54  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_RZDIR_INC
#define CERNLIB_ZEBRA_RZDIR_INC
*
*
* rzdir.inc
*
      PARAMETER (NLPATM=100)
      COMMON /RZDIRN/NLCDIR,NLNDIR,NLPAT
      COMMON /RZDIRC/CHCDIR(NLPATM),CHNDIR(NLPATM),CHPAT(NLPATM)
      CHARACTER*16   CHNDIR,    CHCDIR,    CHPAT
C

#endif
