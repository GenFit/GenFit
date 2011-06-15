*
* $Id: genset.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: genset.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GENSET
*
*      This subroutine closes a SET file, and writes out the
*      end-data assembly and end-SET file blocks.
*      Author :    Nils Hoimyr
*
#include "geant321/gcsetf.inc"
 
*
      WRITE (BLKSTR,10000) N1+1,N1+3
      CALL GJWRIT
* *** the last line
      WRITE (NUNIT1,'(A80)') JVSTR
*
      CLOSE (UNIT = NUNIT2)
      CLOSE (UNIT = NUNIT1)
*
10000 FORMAT('@9998,',I10,'@9999,',I10)
*
      END
