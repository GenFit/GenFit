*
* $Id: cmsres.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: cmsres.inc,v $
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
#ifndef CERNLIB_GEANT321_CMSRES_INC
#define CERNLIB_GEANT321_CMSRES_INC
*
*
* cmsres.inc
*
C***********************************************************************
C     INCLUDE FILE: CMSRES COPY                   CREATED 11/12/86 BY PA
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX5=1133)
      CHARACTER*8 ANR
      COMMON /FKCMSR/PXR(IDMAX5),PYR(IDMAX5),PZR(IDMAX5),HER(IDMAX5),
     *AMR(IDMAX5),ICHR(IDMAX5),IBARR(IDMAX5),            NRER(IDMAX5),
     *NRES
      COMMON /FKCMSC/ANR(IDMAX5)
 

#endif
