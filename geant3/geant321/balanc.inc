*
* $Id: balanc.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: balanc.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_BALANC_INC
#define CERNLIB_GEANT321_BALANC_INC
*
*
* balanc.inc
*
*                                                                     *
*=== balanc ==========================================================*
*                                                                     *
      LOGICAL LRESMP, LEVDIF, LPRDIF
      COMMON /FKBALA/  ETTOT,  PTTOT, PXTTOT, PYTTOT, PZTTOT,
     &                 EINCP,  EINCN, TVGREY, TVGRE0,  TVEUZ,
     &                   EUZ,    PUX,    PUY,    PUZ,  ENUCR, PXNUCR,
     &                PYNUCR, PZNUCR,  EINTR, PXINTR, PYINTR, PZINTR,
     &                  EFRM,  PXFRM,  PYFRM,  PZFRM,   PSEA,
     &                NGREYP, NGREYN,    ICU,    IBU, ICNUCR, IBNUCR,
     &                ICINTR, IBINTR, LRESMP, LEVDIF, LPRDIF
 

#endif
