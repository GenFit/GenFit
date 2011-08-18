      SUBROUTINE UFILES
*
*            To open FFREAD and HBOOK files
*
      OPEN(UNIT=4,FILE='model.dat',STATUS='UNKNOWN')
      END
      SUBROUTINE UGINIT
*
************************************************************************
*                                                                      *
*              To initialise GEANT3 program and read data cards        *
*                                                                      *
************************************************************************
*
#include "gckine.inc"
#include "gcunit.inc"
#include "model.inc"
*
*     -----------------------------------------------------------------
*
*             Open user files
*
      CALL UFILES
*
*             Initialize GEANT
C..geant..
      CALL GINIT
*
*             Prints version number
*
      WRITE(LOUT,1001)
*
*             IKINE  = particle type (default=1=gamma)
*             PKINE(1)=particle energy
*             IKINE and PKINE can be changed with the data card KINE
*
      IKINE=1
      PKINE(1)=10.
*
*             Read data cards with FFREAD
*
C..geant..

      write(lout,1000)
 1000 format(/,' ========> Reading ffread data cards : type <======='
     +,/,'read 4'
     +,/,'your own data cards if any'
     +,/,'stop',/,'      Now waiting for input',/)

      CALL GFFGO
*
*             Initialize GEANT/ZBOOK data structures
*
C..geant..
      CALL GZINIT
*
*             Initialize graphics package
*
      CALL GDINIT
*
*             Geometry and materials description.
*
      CALL UGEOM
*
*             Particle table definition and energy loss initialization.
*
C..geant..
      CALL GPART
C..geant..
      CALL GPHYSI
*
*             Initialize MODE graphics 
*
      CALL VIEWYZ(1)
C..book user histograms:
      CALL UHINIT
*
 1001 FORMAT(/,'  MODE VERSION 1.00 : ',/)
      END

      SUBROUTINE UHINIT
*
************************************************************************
*                                                                      *
*             To book the user's histograms                            *
*                                                                      *
************************************************************************
*
*
#include "gckine.inc"
*
*     ------------------------------------------------------------------
*
      ELOW = 0.
      EHIG = 1.1*PKINE(1)
      CALL HBOOK1(101,'TOT ENERGY IN MODEL$', 100, ELOW, EHIG, 0.)
      CALL HBOOK1(102,'TOT ENERGY IN C6F6$ ', 100, ELOW, EHIG, 0.)
      CALL HBOOK1(103,'TOT ENERGY IN SCIN$ ', 100, ELOW, EHIG, 0.)
*
      END

      SUBROUTINE UGEOM
*
************************************************************************
*                                                                      *
*             Routine to define the geometry of the set-up.            *
*                                                                      *
************************************************************************
*
*
#include "model.inc"
*
      DIMENSION PAR(20)
      DIMENSION AC6F6(2),ZC6F6(2),WC6F6(2)
*
      DATA AC6F6/12.01,19.01/
      DATA ZC6F6/6.,9./
      DATA WC6F6/6.,6./
*
*     -----------------------------------------------------------------
*
*
*             Defines materials
      CALL GSMATE( 1,'AIR$     ',  15.0,7.0,0.0012,30050.0,67500.0,0,0) 
      CALL GSMATE( 2,'PLAST SC$',  6.25,3.4,1.032 ,   43.0,   437.,0,0) 
      CALL GSMATE( 3,'IRON$    ', 55.85,26.,7.8   ,   1.76,   17.1,0,0) 
      CALL GSMATE( 9,'ALUMINIUM$', 26.98,13.,2.7   , 8.9,37.2,0,0)  
      CALL GSMIXT(10,'C6F6$',AC6F6,ZC6F6,1.61,-2,WC6F6)
      CALL GSMATE(13,'LEAD$     ',207.19,82.,11.35 ,0.56,18.5,0,0)  
*
*             Defines tracking media parameters.
      FIELDM =  0.
      IFIELD =  0
      TMAXFD =  0.
      DMAXMS =  0.1
      DEEMAX =  0.02
      EPSIL  =  0.01
      STMIN  =  0.01
*
      CALL GSTMED( 1,'AIR               $'    ,  1 , 0 , IFIELD,    
     *                FIELDM,TMAXFD,0.,0.        , EPSIL, STMIN, 0 , 0 )    
      CALL GSTMED( 2,'PLASTIC SCINTILLAT$'    ,  2 , 0 , IFIELD,    
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )    
      CALL GSTMED( 3,'IRON              $'    ,  3 , 0 , IFIELD,    
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )    
      CALL GSTMED( 9,'ALUM              $'    ,  9 , 0 , IFIELD,    
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )    
      CALL GSTMED(10,'C6F6              $'    , 10 , 0 , IFIELD,
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )
      CALL GSTMED(13,'LEAD              $'    , 13 , 0 , IFIELD,    
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )    
*
*
*             Defines geometry of the set-up
*
*                  Define Mother volume MODE
      PAR(1)=8.5
      PAR(2)=8.5
      PAR(3)=10.
      CALL GSVOLU('MODE', 'BOX ', 1, PAR, 3, IVOLU)
*
*                  Define Iron absorber IROM
      PAR(1)=0.
      PAR(2)=8.
      PAR(3)=4.
      CALL GSVOLU('IROM', 'TUBE', 9, PAR, 3, IVOLU)
      CALL GSPOS('IROM',1,'MODE',0.,0.,-5.8,0,'ONLY')
*
      PAR(1)=0.
      PAR(2)=8.
      PAR(3)=0.4
      CALL GSVOLU('IRON', 'TUBE', 3, PAR, 3, IVOLU)
      DO 10 I=1,8
         Z0=-4.+(I-1)*1.+0.4
         CALL GSPOS('IRON',I,'IROM',0.,0.,Z0,0,'ONLY')
  10  CONTINUE
*      CALL GSORD('IROM',3)
*
*                 Define the 2 Polycones of C6F6
      PAR( 1)=0.
      PAR( 2)=360.
      PAR( 3)=4.
      PAR( 4)=-1.6
      PAR( 5)=0.4
      PAR( 6)=3.8
      PAR( 7)=-1.
      PAR( 8)=0.8
      PAR( 9)=2.8
      PAR(10)=1.
      PAR(11)=0.8
      PAR(12)=2.8
      PAR(13)=2.
      PAR(14)=0.4
      PAR(15)=3.8
      CALL GSVOLU('C6F6','PCON',10, PAR,15, IVOLU)
      CALL GSPOS('C6F6',1,'MODE',-4.,-4.,0.,0,'ONLY')
      CALL GSPOS('C6F6',2,'MODE',-4., 4.,0.,0,'ONLY')
      CALL GSPOS('C6F6',3,'MODE', 4.,-4.,0.,0,'ONLY')
      CALL GSPOS('C6F6',4,'MODE', 4., 4.,0.,0,'ONLY')
*
*                 Define the 2 scintillator plates
      PAR(1)=7.
      PAR(2)=7.
      PAR(3)=0.5
      CALL GSVOLU('SCIN', 'BOX ', 2, PAR, 3, IVOLU)
      CALL GSPOS('SCIN',1,'MODE',0.,0.,3.2,0,'ONLY')
      CALL GSPOS('SCIN',2,'MODE',0.,0.,4.6,0,'ONLY')
*
*                 Define the lead absorber
      PAR(1)=4.
      PAR(2)=7.
      PAR(3)=5.
      PAR(4)=7.
      PAR(5)=1.5
      CALL GSVOLU('LEAD', 'TRD2',13, PAR, 5, IVOLU)
      CALL GSPOS('LEAD',1,'MODE',0.,0.,8.,0,'ONLY')
*
      CALL GSORD('MODE',3)
*
*             Close geometry banks. Mandatory system routine.
*
      CALL GGCLOS
*
      END

      SUBROUTINE GUKINE
*
************************************************************************
*                                                                      *
*             Generates Kinematics for primary tracks                  *
*                                                                      *
************************************************************************
*
#include "gcflag.inc"
#include "gckine.inc"
#include "gconst.inc"
#include "model.inc"
*
      DIMENSION VERTEX(3),PLAB(3)
      DIMENSION RNDM(4)
*
*     -----------------------------------------------------------------
*

      CALL GRNDM(RNDM,4,1)
      XV=8.5*(2.*RNDM(1)-1.)
      YV=8.5*(2.*RNDM(2)-1.)
      THETA=PI*RNDM(3)/10.
      PHI=TWOPI*RNDM(4)
C
      PLAB(1) = PKINE(1)*SIN(THETA)*COS(PHI)
      PLAB(2) = PKINE(1)*SIN(THETA)*SIN(PHI)
      PLAB(3) = PKINE(1)*COS(THETA)
      VERTEX(1)=XV
      VERTEX(2)=YV
      VERTEX(3)=-9.99
C
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
      CALL GSKINE(PLAB,IKINE,NVERT,0,0,NT)
*
*             Reset energy deposited
*
      EC6F6=0.
      ESCIN=0.
      SUMX0=0.
*
*              Kinematic debug (controled by ISWIT(1))
*
      IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
        CALL GPRINT('VERT',0)
        CALL GPRINT('KINE',0)
      ENDIF
*
      END

      SUBROUTINE GUTREV
*
************************************************************************
*                                                                      *
*             User routine to control tracking of one event            *
*                                                                      *
*             Called by GRUN (in GEANT, controls a run of events)      *
*                                                                      *
************************************************************************
*
#include "gcflag.inc"
*
*     -----------------------------------------------------------------
*..geant..(TRAK 110)..
      CALL GTREVE
C..              in turn calls GUTRAK which, if not here,
C..              is in GEANT3. GUTRAK calls GTRACK (steering routine
C..              to track a ptcl).
*
*             Debug and plot tracks.
*
      IF(IDEBUG.EQ.1) THEN
        IF(ISWIT(2).EQ.1) CALL GPRINT('JXYZ', 0)
        IF(ISWIT(3).EQ.1) THEN
C..geant..draw a view..
          CALL GDSHOW(1)
C..geant..draw a track (0 means all tracks)..
          CALL GDXYZ (0)
        ENDIF
      ENDIF
*
      END

      SUBROUTINE GUSTEP
*
************************************************************************
*                                                                      *
*             User routine called at the end of each tracking step     *
*             MEC   is the mechanism origin of the step                *
*             INWVOL is different from 0 when the track has reached    *
*                    a volume boundary                                 *
*             ISTOP is different from 0 if the track has stopped       *
*                                                                      *
************************************************************************
*
#include "gcmate.inc"
#include "gctmed.inc"
#include "gcking.inc"
#include "gcflag.inc"
#include "gctrak.inc"
#include "model.inc"
*
*     -----------------------------------------------------------------
*
*             Accumulate energy deposited in C6F6 and SCINTILLATORS
*
      SUMX0=SUMX0+STEP/RADL
      IF(NUMED.EQ.10)THEN
         EC6F6=EC6F6+DESTEP
      ENDIF
      IF(NUMED.EQ. 2)THEN
         ESCIN=ESCIN+DESTEP
      ENDIF
*
*             Something generated ?
*
      IF(NGKINE.GT.0) THEN
         CALL GSKING(0)
      END IF
*
*             Debug/plot event
      IF(IDEBUG.EQ.1) THEN
        IF((ISWIT(2).EQ.1).OR.(ISWIT(3).EQ.1)) CALL GSXYZ
        IF (ISWIT(2).EQ.2) CALL GPCXYZ
        IF (ISWIT(2).EQ.3) CALL GDCXYZ
      ENDIF
      END

      SUBROUTINE GUOUT
*
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine called at the end of each event.            *
C.    *                                                                *
C.    ******************************************************************
C.
C.
#include "gcunit.inc"
#include "gcflag.inc"  
#include "model.inc"
*
      DATA ID101,ID102,ID103,ID104/4*0/
C.
C.    ------------------------------------------------------------------
C.
      ETOT=EC6F6+ESCIN
      IF(ETOT.GT.0.)THEN
         CALL HFF1(101,ID101,ETOT,1.)
         IF(EC6F6.GT.0.)CALL HFF1(102,ID102,EC6F6,1.)
         IF(ESCIN.GT.0.)CALL HFF1(103,ID103,ESCIN,1.)
      ENDIF
C
      IF(IDEBUG.NE.0)THEN
         IF(ISWIT(8).NE.0)THEN
            PRINT 1000,IEVENT,XV,YV,THETA,PHI,SUMX0
 1000       FORMAT(' IEVENT=',I6,' XV=',F8.2,' YV=',F8.2,
     +             ' THETA=',F8.2,' PHI=',F8.2,' SUMX0=',F10.4)
         ENDIF
         IF(ISWIT(6).NE.0)THEN
            WRITE(LOUT,701)IEVENT,ETOT,EC6F6,ESCIN
 701        FORMAT(' - SUMMARY OF EVENT NR ',I6,
     +      ' ENERGY DEPOSIT IN MODE :',E11.3,' GEV, ',
     +      ' IN C6F6 :',E11.3,
     +      ' IN SCIN :',E11.3)
*
         ENDIF
      ENDIF
  100 RETURN
*
      END
C
      SUBROUTINE UGLAST
*
************************************************************************
*                                                                      *
*            Termination routine to print histograms and statistics    *
*                                                                      *
************************************************************************
*
*     -----------------------------------------------------------------
*
      CALL GLAST
*
*             Print HBOOK histograms
*
      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('HBOOK')
      CALL HPRINT(0)
*
      END
      SUBROUTINE VIEWYZ (IVIEW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Draw full set up in 'view bank' mode.                      *
C.    *                                                                *
C.    ******************************************************************
C.
C.
C.    ------------------------------------------------------------------
C.
C             Create bank for view XY.
C
      CALL GDOPEN(IVIEW)
C
      CALL GDHEAD(110110, 'VIEW XY$' ,0.5)
C
      CALL GDRAWC ('MODE', 1,0., 10.,10.,1.,1.)
C
      END
