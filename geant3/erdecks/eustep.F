*
* $Id: eustep.F 227 2008-06-09 12:09:52Z brun $
*
* $Log: eustep.F,v $
* Revision 1.2  2007/05/18 08:44:16  brun
* A major update of GEANTE by Andrea Fontana and Alberto Rotondi
*
*
* 1) update of the Coulomb multiple scattering parametrization;
* 2) update of the straggling of energy loss for thin materials;
* 3) new options to extrapolate the track parameters to the point
*    of closest approach to a point or to a wire (straight line).
*
* Details on the physical motivation behind this work can be found
* in our report for the Panda Collaboration, available at:
*
* http://www.pv.infn.it/~fontana/tracking.pdf
*
* Feel free to contact us for questions and discussions about these
* topics by using the following email addresses:
*
* alberto.rotondi@pv.infn.it
* andrea.fontana@pv.infn.it
*
* ---
*
* List of changes in the fortran and C++ routines of the geant3
* VMC directory:
*
* - gcmore.inc
*   gtmore.inc
*   geant3LinkDef.h
*   gcomad.F
*
*  Added a new common that contains all the new variables:
*       COMMON/GCMORE/GCALPHA,ICLOSE,PFINAL(3),DSTRT,WIRE1(3),WIRE2(3),
*      +              P1(3),P2(3),P3(3),CLENG(3)
*
*      input to ERLAND:
*       GCALPHA: energy cut parameter for energy loss fluctuations
*
*      input to EUSTEP:
*       ICLOSE: = 1 the use of the common is enabled for the closest
*                   approach to a point PFINAL(3)
*               = 2 the use of the common is enabled for the closest
*                   approach to a wire of extremes WIRE1(3) and WIRE2(3)
*               = 0 the common is empty and disabled
*       PFINAL(3): assigned point
*       DSTRT: assigned distance between initial point in ERTRAK
*              and PFINAL along straight line (currently noy used)
*       WIRE1(3): first point of a wire
*       WIRE2(3): second point of a wire
*
*      output from EUSTEP:
*       P1(3): point previous to the point of closest approach to
*              PFINAL() or wire
*       P2(3): point of closest approach to PFINAL() or wire
*       P3(3): point next to the point of closest approach to
*              PFINAL() or wire
*       CLENG(3): track length to the previous 3 points
*
*       Important note: the calculated points of closest approach are
*       depending on the GEANE steps. For calculating the true point
*       of closest approach the last 3 points of the extrapolation, i.e.
*       the previous to closest, the closest and the next to closest are
*       returned to the user. Different algorithms can be implemented, but
*       we decided to leave this to the users in the C++ interface to GEANE.
*
* - ermcsc.F
*  new expression for the variance of the Coulomb multiple scattering
*  according to Fruhwirth and Regler, NIM A 456 (2001) 369
*
* - ertrch.F
*  added DESTEP in the calling string of ERLAND for calculation with
*  Urban model. Added and saved previous step PRSTEP.
*
* - erland.F
*  added new calculation for sigma of straggling in energy loss
*  to include in Geane the Urban/Landau approximation, as explained
*  in the Geant manual and related papers.
*  The model parametrization can be controlled with a user variable (GCALPHA)
*  in the new GCMORE common block: 1 is for old gaussian model valid
*  for dense materials, other values (see the report) are for gaseous
*  materials.
*
* - eustep.F
*  added the calculation to the distance of closest approach to a point
*  or to a wire.
*
* - TGeant3.h
* - TGeant3.cxx
*  added the possibility to define user cuts (already present in the gccuts
*  struct but not in the TGeant3::SetCUTS method) and to define the new
*  variables of the GCMORE common with two new methods SetECut() and
*  SetClose().
*  Added new method InitGEANE() to initialize GEANE to the old behaviour
*  (default) for backward compatibility. Only the multiple scattering has
*  been updated to a more correct formula.
*  Corrected a typo in the call to the routine Trscsd().
*
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:36  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE EUSTEP
* modified by A.Rotondi & A.Fontana 
* may 2007
* (to be used to find the closest approach to a point
*  and to a wire)
*
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconst.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcmore.inc"
#include "geant321/gcjump.inc"
*
      INTEGER*4 ICC,CFLAG
      SAVE ICC, ICONT
      REAL*4 DIST2,PRDIST2
  
      SAVE PRDIST2
      REAL*4 D2X,D2Y,D2Z,AMODD
      DATA ICC/0/, ICONT/0/
*
*     clear ICC when tracking starts
      IF(SLENG .EQ. 0.) ICC = 0
     
      CFLAG=ICLOSE
      IF (IDEBUG*ISWIT(2).NE.0) CALL ERXYZC
      
      IF(CFLAG.EQ.1) THEN
*      distance between the track point and the point        
         IF(ICC.EQ.0) PRDIST2=BIG
         DIST2 = (VECT(1)-PFINAL(1))**2+(VECT(2)-PFINAL(2))**2
     +        + (VECT(3)-PFINAL(3))**2
         
         IF((SQRT(DIST2)-SQRT(PRDIST2)).LT.1.E-3) THEN
            PRDIST2 = DIST2
            ICC = 1
            ICONT = 1
            CLENG(1) = CLENG(2)
            CLENG(2) = CLENG(3)   
            CLENG(3) = SLENG
            CALL UCOPY(P2,P1,3)
            CALL UCOPY(P3,P2,3)
            CALL UCOPY(VECT,P3,3)
         ELSE
*      store the first point of increasing distance
            IF(ICONT .EQ. 1) THEN
               CLENG(1) = CLENG(2)
               CLENG(2) = CLENG(3)   
               CLENG(3) = SLENG
               CALL UCOPY(P2,P1,3)
               CALL UCOPY(P3,P2,3)
               CALL UCOPY(VECT,P3,3)
               ICONT = 0
            ENDIF
         ENDIF
         
      ELSE IF(CFLAG.EQ.2) THEN
         
         IF(ICC.EQ.0) PRDIST2=BIG
         
         D2X = (WIRE2(2)-WIRE1(2))*(WIRE1(3)-VECT(3))
     +        - (WIRE2(3)-WIRE1(3))*(WIRE1(2)-VECT(2)) 
         D2Y = (WIRE2(3)-WIRE1(3))*(WIRE1(1)-VECT(1))
     +        - (WIRE2(1)-WIRE1(1))*(WIRE1(3)-VECT(3)) 
         D2Z = (WIRE2(1)-WIRE1(1))*(WIRE1(2)-VECT(2))
     +        - (WIRE2(2)-WIRE1(2))*(WIRE1(1)-VECT(1)) 
         AMODD = (WIRE2(1)-WIRE1(1))**2 + (WIRE2(2)-WIRE1(2))**2 
     +        + (WIRE2(3)-WIRE1(3))**2   
         DIST2 = (D2X**2 + D2Y**2 + D2Z**2)/AMODD
*      distance between the track point and the wire   
         IF((SQRT(DIST2)-SQRT(PRDIST2)).LT.1.E-3) THEN
            PRDIST2 = DIST2
            ICC=1
            ICONT = 1	    
            CLENG(1) = CLENG(2)
            CLENG(2) = CLENG(3)   
            CLENG(3) = SLENG
            CALL UCOPY(P2,P1,3)
            CALL UCOPY(P3,P2,3)
            CALL UCOPY(VECT,P3,3)
         ELSE
*      store the first point of increasing distance
            IF(ICONT .EQ. 1) THEN
               CLENG(1) = CLENG(2)
               CLENG(2) = CLENG(3)   
               CLENG(3) = SLENG
               CALL UCOPY(P2,P1,3)
               CALL UCOPY(P3,P2,3)
               CALL UCOPY(VECT,P3,3)
               ICONT = 0
            ENDIF
         ENDIF
      ENDIF
      CALL GUSTEP
*     
      END
