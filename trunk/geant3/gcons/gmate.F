*
* $Id: gmate.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: gmate.F,v $
* Revision 1.1.1.1  2002/07/24 15:56:24  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:37  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:16  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE G3MATE
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Define standard GEANT materials                               *
C.    *                                                                *
C.    *  All data EXCEPT for nuclear absorption lengths taken from :   *
C.    *                        M. Aguilar-Benitez et al,               *
C.    *                        Rev. of Particle Properties,            *
C.    *                        Rev. Mod. Phys. 56(1984)                *
C.    *                                                                *
C.    *    ==>Called by : <USER>, UGINIT                               *
C.    *       Author    G.Patrick  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION UBUF(1)
      DATA UBUF/0./
C
      CALL G3SMATE( 2,'DEUTERIUM$',  2.01, 1.,0.162 ,757.,342.,UBUF,0)
      CALL G3SMATE( 3,'HELIUM$   ',  4.  , 2.,0.125 ,755.,478.,UBUF,0)
      CALL G3SMATE( 4,'LITHIUM$  ',  6.94, 3.,0.534 ,155.,120.6,UBUF,0)
      CALL G3SMATE( 5,'BERILLIUM$',  9.01, 4.,1.848 ,35.3,36.7,UBUF,0)
      CALL G3SMATE( 6,'CARBON$   ', 12.01, 6.,2.265 ,18.8,49.9,UBUF,0)
      CALL G3SMATE( 7,'NITROGEN$ ', 14.01, 7.,0.808 ,44.5,99.4,UBUF,0)
      CALL G3SMATE( 8,'NEON$     ', 20.18,10.,1.207 , 24.,74.9,UBUF,0)
      CALL G3SMATE( 9,'ALUMINIUM$', 26.98,13.,2.7   , 8.9,37.2,UBUF,0)
      CALL G3SMATE(10,'IRON$     ', 55.85,26.,7.87  ,1.76,17.1,UBUF,0)
      CALL G3SMATE(11,'COPPER$   ', 63.54,29.,8.96  ,1.43,14.8,UBUF,0)
      CALL G3SMATE(12,'TUNGSTEN$ ',183.85,74.,19.3  ,0.35,10.3,UBUF,0)
      CALL G3SMATE(13,'LEAD$     ',207.19,82.,11.35 ,0.56,18.5,UBUF,0)
      CALL G3SMATE(14,'URANIUM$  ',238.03,92.,18.95 ,0.32,12. ,UBUF,0)
      CALL G3SMATE(15,'AIR$      ',14.61,7.3,0.001205,30423.24,67500.
     +                           ,UBUF,0)
      CALL G3SMATE(16,'VACUUM$ ',1.E-16,1.E-16,1.E-16,1.E+16,1.E+16
     +                           ,UBUF,0)
C
      END
