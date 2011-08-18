*
* $Id: gufld.F 233 2009-01-21 14:47:44Z brun $
*
* $Log: gufld.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:42  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GUFLD(X,F)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute the magnetic field F             *
C.    *       at space point X                                         *
C.    *                                                                *
C.    *    ==>Called by : G3HELIX,G3RKUTA                              *
C.    *                                                                *
C.    *    Please note that starting from svn revision 232, the        *
C.    *    argument of this routine have been moved from single to     *
C.    *    double precision. The corresponding changes have been       *
C.    *    done to the calling routines. This change is backward       *
C.    *    incompatible and to continue working with single precision  *
C.    *    arguments, the compilation flag -DSINGLEFIELD should be     *
C.    *    added to the compilation, both FORTRAN and C++.             *
C.    *                                                                *
C.    *    F.Carminati, January 20, 2009                               *
C.    *                                                                *
C.    ******************************************************************
C.
#ifdef SINGLEFIELD
      DIMENSION X(3),F(3)
#else
      DOUBLE PRECISION X(3),F(3)
#endif
C.
C.    ------------------------------------------------------------------
C.
      F(1)=0.
      F(2)=0.
      F(3)=0.
C
      END
