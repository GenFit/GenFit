*
* $Id: imp128.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: imp128.inc,v $
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
* imp128.inc
*
#if (!defined(CERNLIB_F90))&&(defined(CERNLIB_NOQUAD)||!defined(CERNLIB_DOUBLE))
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#elif (defined(CERNLIB_F90))&&(!defined(CERNLIB_DOUBLE))
      IMPLICIT REAL(2) (A-H,O-Z)

#elif (defined(CERNLIB_F90))&&(defined(CERNLIB_DOUBLE))
      IMPLICIT REAL(4) (A-H,O-Z)

#elif 1
      IMPLICIT REAL*16 (A-H,O-Z)

#endif
