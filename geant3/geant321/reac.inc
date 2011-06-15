*
* $Id: reac.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: reac.inc,v $
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
#ifndef CERNLIB_GEANT321_REAC_INC
#define CERNLIB_GEANT321_REAC_INC
*
*
* reac.inc
*
*
*=== reac =============================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 10 december 1991  by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 10-dec-91     by    Alfredo Ferrari               *
*                                                                      *
*     This is the original common reac of Hadrin                       *
*                                                                      *
*     Included in the following routines:                              *
*                                                                      *
*                CALUMV                                                *
*                DATESH                                                *
*                HADRIV                                                *
*                RCHANV                                                *
*                BLKDT3                                                *
*                HADRIN                                                *
*                IEFUN                                                 *
*                SIGINT                                                *
*                CALUMO                                                *
*                RCHANW                                                *
*                                                                      *
*----------------------------------------------------------------------*
*
      COMMON /FKREAC/ UMO   ( 296), PLABF ( 296), SIIN  ( 296),
     +              WK    (5184), NRK  (2,268), NURE  (30,2)
 

#endif
