*
* $Id: imp64.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:36  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
#if !defined(CERNLIB_DOUBLE)
      IMPLICIT REAL (A-H,O-Z)

#elif (defined(CERNLIB_F90))&&(defined(CERNLIB_DOUBLE))
      IMPLICIT REAL(2) (A-H,O-Z)

#elif 1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#endif
