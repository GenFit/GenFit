*
* $Id: part2.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: part2.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:42  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_PART2_INC
#define CERNLIB_GEANT321_PART2_INC
*
*
* part2.inc
*
*
*=== part2 ============================================================*
*
*----------------------------------------------------------------------*
*     Include file: part2 copy        Revised on 20-7-90 by A. Ferrari *
*     Note: see also part copy and part3 copy                          *
*     Changes: none                                                    *
*     Included in the following subroutines or functions: not updated  *
*                                                                      *
*     Description of the common block(s) and variable(s)               *
*                                                                      *
*         Kptoip = conversion from part to paprop numbering            *
*         Iptokp = conversion from paprop to part numbering            *
*                                                                      *
*----------------------------------------------------------------------*
*
      PARAMETER (IDMAX8=183)
      CHARACTER*8  ANAME
      COMMON / FKPART / AAM  (IDMAX8), GA     (IDMAX8), TAU    (IDMAX8),
     &                AAMDSC (IDMAX8), IICH   (IDMAX8), IIBAR  (IDMAX8),
     &                K1     (IDMAX8), K2     (IDMAX8), KPTOIP (IDMAX8),
     &                IPTOKP (39)
      COMMON / FKCHPA / ANAME (IDMAX8)
 

#endif
