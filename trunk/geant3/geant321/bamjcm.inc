*
* $Id: bamjcm.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: bamjcm.inc,v $
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
#ifndef CERNLIB_GEANT321_BAMJCM_INC
#define CERNLIB_GEANT321_BAMJCM_INC
*
*
* bamjcm.inc
*
*
*=== bamjcm ===========================================================*
*
 
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Bamjcm by A. Ferrari                                             *
*                           included in:                               *
*                                        bamjet                        *
*     When changing kmxjcm dimension look also at verein!!!!!          *
*----------------------------------------------------------------------*
*                                                                      *
      PARAMETER (KMXJCM = 100)
      COMMON /FKBAMJ/  RPX(KMXJCM),  RPY(KMXJCM),   RE(KMXJCM),
     &                KFR1(KMXJCM), KFR2(KMXJCM),   IV(KMXJCM)
 
 

#endif
