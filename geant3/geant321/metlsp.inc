*
* $Id: metlsp.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: metlsp.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:40  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_METLSP_INC
#define CERNLIB_GEANT321_METLSP_INC
*
*
* metlsp.inc
*
C***********************************************************************
C     INCLUDE FILE: METLSP COPY                   CREATED 11/12/86 BY PA
C     CHANGES: NONE
C     INCLUDED IN THE FOLLOWING SUBROUTINES OR FUNCTIONS: NOT UPDATED
C
C     DESCRIPTION OF THE COMMON BLOCK(S) AND VARIABLE(S)
C
C***********************************************************************
      PARAMETER (IDMAX2=1133)
      COMMON/FKMETL/ IS,ITS(IDMAX2),CXS(IDMAX2),CYS(IDMAX2),CZS(IDMAX2),
     *ELS(IDMAX2),PLS(IDMAX2)
 

#endif
