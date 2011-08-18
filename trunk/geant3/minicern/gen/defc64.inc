*
* $Id: defc64.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: defc64.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:46  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:35  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* defc64.inc
*
#if (defined(CERNLIB_DOUBLE))&&(defined(CERNLIB_F90))
      COMPLEX(2)

#elif !defined(CERNLIB_DOUBLE)
      COMPLEX

#elif 1
      COMPLEX*16

#endif
