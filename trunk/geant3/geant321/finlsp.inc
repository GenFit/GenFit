*
* $Id: finlsp.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: finlsp.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:38  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_FINLSP_INC
#define CERNLIB_GEANT321_FINLSP_INC
*
*
* finlsp.inc
*
C***********************************************************************
C     INCLUDE FILE: FINLSP COPY                   CREATED 11/12/86 BY PA
C     NOTE: SEE ALSO FILE FINLSP2 COPY
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX7=100)
      COMMON/FKFINL/ CXR(IDMAX7),CYR(IDMAX7),CZR(IDMAX7),
     *ELR(IDMAX7),PLR(IDMAX7),IR,ITR(IDMAX7)
 

#endif
