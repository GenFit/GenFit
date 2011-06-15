*
* $Id: s_csdim.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: s_csdim.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:43  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_S_CSDIM_INC
#define CERNLIB_GEANT321_S_CSDIM_INC
*
*
* s_csdim.inc
*
C --- DIMENSION STATEMENTS FOR CROSS SECTION DATA ---
      DIMENSION PLAB(41),CSEL(35,41),CSIN(35,41),CSPIEL(3,41),
     $          CSPIIN(3,41),CSPNEL(3,41),CSPNIN(3,41),
     $          ELAB(17),CNLWAT(15),CNLWEL(15,17),CNLWIN(15,17),
     $          CSCAP(100),EKFISS(21),CSFISS(4,21)
C

#endif
