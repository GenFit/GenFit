      SUBROUTINE GUKINE
*
************************************************************************
*                                                                      *
*             GEANT3 user routine to generate Kinematics               *
*               for primary tracks                                     *
*                                                                      *
************************************************************************
*
#include "gckine.inc"  
#include "gcflag.inc"  
#include "gconsp.inc"  
*KEEP,URGEOM.
*
      COMMON/URGEOM/UTHICK,CSHIEL,ESHIEL,BTHICK,AIRGP1,CTHICK,AIRGP2,
     *              AIRGP3,EPCHAM,HMODL1,HSECT1,HMODL2,HSECT2,HMODL3,
     *              HSECT3,HLCALO,WGS1,WGS2,IGAST
C
*KEEP,URMIPS.
      COMMON/URMIPS/XMIP,NTXMIP,TXMIP,TXMIP2,CALIB1,CALIB2,XKPL(60,8),
     *              XPHMIP(20,34),IDSTW
C
*KEND.
      DIMENSION VERTEX(6),PLAB(3),RNDM(2)
*
*     -----------------------------------------------------------------
*
      XMIP=0.
      CALL VZERO(VERTEX,6)
      VERTEX(1)=0.04
      VERTEX(2)=0.04
      VERTEX(3)=-HLCALO+.01
      IF(IKINE.GT.100)THEN
         IK=IKINE-100
         THETA=PKINE(2)*DEGRAD
         PHI=PKINE(3)*DEGRAD
         if(pkine(4).ne.0.)vertex(1)=pkine(4)
         if(pkine(5).ne.0.)vertex(2)=pkine(5)
         if(pkine(6).ne.0.)vertex(3)=pkine(6)
      ELSEIF(IKINE.eq.1)then
         inumb=pkine(1)
         in=0
  10     read(3,1000,end=80)itype,ntr,vertex(1),vertex(2),vertex(3)
         if(itype.ne.0)go to 10
         in=in+1
         if(inumb.ne.0)then
            if(in.ne.inumb)go to 10
         endif
         call gsvert(vertex,0,0,0,0,nvert)
         do 20 i=1,ntr
            read(3,1000,end=80)itype,ik,plab
            if(itype.ne.i)go to 10
            call gskine(plab,ik,nvert,0,0,nt)
  20     continue
         go to 90
  80     rewind(3)
         go to 10
      ELSEIF(IKINE.eq.2)then
         inumb=pkine(1)
         in=0
  81     continue
         read(23,1000,end=89)itype,nv,aa,bb,cc
         if(itype.ne.-1)go to 81
         in=in+1
         if(inumb.ne.0)then
            if(in.ne.inumb)go to 81
         endif
         do 85 iv=1,nv
  82        read(23,1000,end=89)itype,ntr,vertex(1),vertex(2),vertex(3)
            if(itype.ne.0)go to 82
            call gsvert(vertex,0,0,0,0,nvert)
            do 84 i=1,ntr
               read(23,1000)itype,ik,plab
               call gskine(plab,ik,nvert,0,0,nt)
  84        continue
  85     continue
         go to 90
  89     rewind(23)
         go to 81
      ELSE
         IK=IKINE
         CALL GRNDM(RNDM,2)
         THETA=PI*RNDM(1)
         PHI=TWOPI*RNDM(2)
         VERTEX(1)=0
         VERTEX(2)=0
         VERTEX(3)=PKINE(3)
      ENDIF
C
      PLAB(1) = PKINE(1)*SIN(THETA)*COS(PHI)
      PLAB(2) = PKINE(1)*SIN(THETA)*SIN(PHI)
      PLAB(3) = PKINE(1)*COS(THETA)
C
      CALL GSVERT(VERTEX,0,0,0,0,NVERT)
      CALL GSKINE(PLAB,IK,NVERT,0,0,NT)
*
*              Kinematic debug (controled by ISWIT(1))
*
  90  IF(IDEBUG.EQ.1.AND.ISWIT(1).EQ.1) THEN
        CALL GPRINT('VERT',0)
        CALL GPRINT('KINE',0)
      ENDIF
*
 1000 format(2i5,3(2x,e14.7))
      END


      SUBROUTINE GUOUT
*
************************************************************************
*                                                                      *
*             GEANT3 user routine called at the end of each event.     *
*                                                                      *
************************************************************************
*
*
#include "gcunit.inc"  
#include "gcflag.inc"  
*KEEP,URMIPS.
      COMMON/URMIPS/XMIP,NTXMIP,TXMIP,TXMIP2,CALIB1,CALIB2,XKPL(60,8),
     *              XPHMIP(20,34),IDSTW
C
*KEND.
      SAVE N10
      DATA N10/0/
*
*     ------------------------------------------------------------------
*
      NTXMIP=NTXMIP+1
      TXMIP=TXMIP+XMIP
      TXMIP2=TXMIP2+XMIP*XMIP
      CALL HFF1(10,N10,XMIP,1.)
*
      IF(ISWIT(6).NE.0)THEN
         IS6=ISWIT(6)
         IF(MOD(IEVENT,IS6).EQ.0)CALL PXMIPS
      ENDIF
*
      END


      SUBROUTINE GUPHAD
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute Hadron. inter. probabilities     *
C.    *                                                                *
C.    *    ==>Called by : GTHADR,GTNEUT                                *
C.    *                                                                *
C.    ******************************************************************
C.
C.
C.    ------------------------------------------------------------------
C.
#include "gcphys.inc"  
*KEND.
C      

      IF (IHADR.NE.4) THEN
         CALL GPGHEI
      ELSE
         CALL FLDIST
      ENDIF
      END

      SUBROUTINE GUHADR
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to generate one hadronic interaction        *
C.    *                                                                *
C.    *    ==>Called by : GTHADR,GTNEUT                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gcphys.inc"  
C.
C.    ------------------------------------------------------------------
C.
C
C          GHEISHA only if IHADR<3 (default)
C
      IF (IHADR.NE.4) THEN
         CALL GHEISH
      ELSE
         CALL FLUFIN
      ENDIF
      END


      SUBROUTINE GUSTEP
*
************************************************************************
*                                                                      *
*        GEANT3 user routine called at the end of each tracking step   *
*                                                                      *
************************************************************************
*
#include "gcbank.inc"  
#include "gctmed.inc"  
#include "gckine.inc"  
#include "gcking.inc"  
#include "gcflag.inc"  
#include "gctrak.inc"  
#include "gcvolu.inc"  
#include "gcscan.inc"  
*KEEP,URMIPS.
      COMMON/URMIPS/XMIP,NTXMIP,TXMIP,TXMIP2,CALIB1,CALIB2,XKPL(60,8),
     *              XPHMIP(20,34),IDSTW
C
*KEND.
      DIMENSION X(13)
      SAVE KMOD2,KMOD3,N500,N5,N6,XMIPMA,TOMAX,IPLOLD,PLMIP,SSTEP,FIRST
      CHARACTER*4 CMOD2, CMOD3
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CMOD2, CMOD3 / 'MOD2', 'MOD3' /
      DATA N500/0/
      DATA N5/0/
      DATA N6/0/
      DATA XMIPMA/70./
      DATA TOMAX/2.E-7/
*
*     -----------------------------------------------------------------
*
      IF(FIRST) THEN
         CALL UCTOH(CMOD2, KMOD2, 4, 4)
         CALL UCTOH(CMOD3, KMOD3, 4, 4)
         FIRST = .FALSE.
      ENDIF
      IF(SLENG.LE.0.)THEN
         IPLOLD=0
         PLMIP=0.
         SSTEP=0.
      ENDIF
*             Something generated ?
      IF(NGKINE.GT.0) THEN
        DO 5 I=1,NGKINE
          ITYPA  = GKIN(5,I)
          UPWGHT = ITYPA
          IF(ITYPA.NE.4)  CALL GSKING(I)
   5    CONTINUE
      ENDIF
*
*             Are we in the gas ?
*
      IF(NUMED.EQ.5)THEN
*
*             Evaluate average energy loss for muons
*
         IF(ITRTYP.EQ.5)THEN
            IF(DESTEP.GT.0..AND.GETOT.GT.0.3)THEN
               IF(STEP.GT.0.1)THEN
                  DEDX=DESTEP/STEP
                  CALL HFF1(500,N500,DEDX,1.)
               ENDIF
            ENDIF
         ENDIF
*
*             Compute number of MIPS
*
         IPL=NUMBER(3)
*         if (ipl.eq.0)ipl=1
         IF(NAMES(3).EQ.KMOD2)IPL=IPL+64
         IF(NAMES(3).EQ.KMOD3)IPL=IPL+99
         IF(DESTEP.GT.0.)THEN
            IF(IPL.LE.64)THEN
               DXMIP=DESTEP/CALIB1
            ELSE
               DXMIP=DESTEP/CALIB2
            ENDIF
            IF(IPL.EQ.IPLOLD)THEN
               PLMIP=PLMIP+DXMIP
               SSTEP=SSTEP+STEP
            ELSE
               PLMIP=DXMIP
               SSTEP=STEP
               IPLOLD=IPL
            ENDIF
            CALL HFF1(5,N5,DXMIP,1.)
            CALL HFF1(6,N6,DXMIP,1.)
            IF(PLMIP.LT.XMIPMA.AND.TOFG.LT.TOMAX)THEN
               XMIP=XMIP+DXMIP
            ENDIF
         ENDIF
*
*             Count number of particles crossing planes
*
         IF(INWVOL.EQ.1)THEN
            IF(MOD(IPL,2).EQ.0)THEN
               KPL=IPL/2
               XKPL(KPL,ITRTYP)=XKPL(KPL,ITRTYP)+1.
            ENDIF
         ENDIF
      ENDIF
*
*             Debug/plot event
*
      if(idebug.ne.0.and.inwvol.ne.2)then
            X(1)=vect(1)
            X(2)=vect(2)
            X(3)=vect(3)
            X(4)=gekin
            X(5)=tofg
            X(6)=istop+1000*ipart
            if(sleng.le.0.)x(6)=x(6)+1000000.
            if(iswit(7).ne.7)then
               CALL HFN(1,X)
            else
               if(sleng.gt.0..and.inwvol.eq.1)call hfn(1,x)
            endif
      endif
      CALL GDEBUG
*     IF(ISTOP.NE.0) THEN
*        IF(IPART.EQ.13) THEN
*          CALL GDTRAK('R')
*        ELSE
*          CALL GDTRAK('DR')
*        ENDIF
*     ENDIF
*
      END


      SUBROUTINE GUTREV
*
************************************************************************
*                                                                      *
*             GEANT3 user routine to control tracking of one event     *
*                                                                      *
*             Called by GRUN                                           *
*                                                                      *
************************************************************************
*
#include "gcflag.inc"  
*
*     -----------------------------------------------------------------
*
      CALL GTREVE
*
*             Debug and plot tracks.
*
      IF(IDEBUG.EQ.1) THEN
        IF(ISWIT(2).EQ.1) CALL GPRINT('JXYZ', 0)
      ENDIF
*
      END


      SUBROUTINE PXMIPS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Print MIPS                                                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gcflag.inc"  
*KEEP,URMIPS.
      COMMON/URMIPS/XMIP,NTXMIP,TXMIP,TXMIP2,CALIB1,CALIB2,XKPL(60,8),
     *              XPHMIP(20,34),IDSTW
C
*KEND.
C
C              Compute average number of MIPS
C
      IF(NTXMIP.GT.0)THEN
         TXM=TXMIP/NTXMIP
         TXMRMS=SQRT(MAX(TXMIP2/NTXMIP - TXM*TXM,0.))
         PRINT 1000,IEVENT,XMIP,TXM,TXMRMS
 1000    FORMAT(' =====> IEVENT=',I5,' MIPS =',F7.1,
     +          '  TOTAL MIPS =',F7.1,' +-',F7.1)
      ENDIF
C
      END


      SUBROUTINE UFILES
*
*            To open FFREAD and HBOOK files
*
      CHARACTER*(*) FILNAM, FSTAT
      PARAMETER (FILNAM='gexam4.dat')
*
      PARAMETER (FSTAT='OLD')
      common/quest/iquest(100)
*
c      open(unit=3,file='opalev.dat',status='old')
c      open(unit=23,file='atlas.dat',status='old')
      OPEN(UNIT=4,FILE=FILNAM,STATUS=FSTAT,
     +     FORM='FORMATTED')
*
*             Open a HBOOK direct access file
*
      iquest(10)=64000
      CALL HROPEN(34,'HBOOK','gexam4.hist','NQ',1024,ISTAT)
      END


      SUBROUTINE UGEOM
*
************************************************************************
*                                                                      *
*             Routine to define the geometry of the set-up.            *
*                                                                      *
************************************************************************
*
#include "gconsp.inc"  
*KEEP,URGEOM.
*
      COMMON/URGEOM/UTHICK,CSHIEL,ESHIEL,BTHICK,AIRGP1,CTHICK,AIRGP2,
     *              AIRGP3,EPCHAM,HMODL1,HSECT1,HMODL2,HSECT2,HMODL3,
     *              HSECT3,HLCALO,WGS1,WGS2,IGAST
C
*KEND.
*
      DIMENSION PAR(10),ABRASS(2),ZBRASS(2),WBRASS(2)
      DIMENSION AISOBU(3),ZISOBU(3),WGAS(3),ACO2(3),ZCO2(3)
      DIMENSION AURAN(2),ZURAN(2),WURAN(2)
      DATA AISOBU/12.01,1.01,39.95/
      DATA ZISOBU/6.,1.,18./
      DATA ACO2/12.01,16.01,39.95/
      DATA ZCO2/6.,8.,18./
      DATA ABRASS/63.54,65.37/
      DATA ZBRASS/29.,30./
      DATA WBRASS/.7,.3/
C
C             Uranium mixture
C
      DATA AURAN/235.,238./
      DATA ZURAN/92.,92./
      DATA WURAN/0.004,0.996/
*
*     -----------------------------------------------------------------
*
*
*             Defines materials
      CALL GSMATE( 1,'AIR$     ',  15.0,7.0,0.0012,30050.0,67500.0,0,0)
      CALL GSMATE( 2,'COPPER$  ', 63.54,29.,8.960 ,   1.43,   14.8,0,0)
      CALL GSMIXT( 3,'URANIUM$  ',AURAN,ZURAN,18.95 ,2,WURAN)
      CALL GSMATE( 4,'CARBON$  ', 12.01, 6., 2.265,  18.80,   49.9,0,0)
      CALL GSMIXT( 6,'BRASS$    ',ABRASS,ZBRASS,8.560   , 2,WBRASS)
*
      IF(IGAST.EQ.1)THEN
*
*             Argon/Isobuthane mixture (60% Ar and 40% Isobuthane)
*             First define Isobuthane compound and relative weights
*
         DISO   =0.00267
         DENS1  =0.002136
         WGAS(1)=4.
         WGAS(2)=10.
         CALL GSMIXT(15,'ISOBUTHAN$',AISOBU,ZISOBU,0.40*DISO,-2,WGAS)
         WGAS(1)=0.40*WGAS(1)
         WGAS(2)=0.40*WGAS(2)
         WGAS(3)=0.60
         CALL GSMIXT( 5,'ARG/ISOBU$',AISOBU,ZISOBU,DENS1, 3,WGAS)
      ELSE
*
*             Argon/CO2 mixture (10% Ar and 90% CO2)
*             First define CO2 compound and relative weights
*
         DENS2  =0.0019573
         DCO2   =0.001977
         WGAS(1)=1.
         WGAS(2)=2.
         CALL GSMIXT(15,'CO2$',ACO2,ZCO2,0.90*DCO2,-2,WGAS)
         WGAS(1)=0.90*WGAS(1)
         WGAS(2)=0.90*WGAS(2)
         WGAS(3)=0.10
         CALL GSMIXT( 5,'ARG/CO2$',ACO2,ZCO2,DENS2, 3,WGAS)
      ENDIF
*
*             Defines tracking media parameters.
      STEMAX =  BIG
      DEEMAX =  0.30
      EPSIL  =  0.005
      STMIN  =  0.1
*
      CALL GSTMED( 1, 'AIR$     ', 1, 0, 0, 0., 0., STEMAX,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
*
      CALL GSTMED( 2, 'COPPER$  ', 2, 0, 0, 0., 0., STEMAX,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
*
      CALL GSTMED( 3, 'URANIUM$ ', 3, 0, 0, 0., 0., STEMAX,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
*
      CALL GSTMED( 4,'CARBON$   ', 4, 0, 0, 0., 0., STEMAX,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
*
      CALL GSTMED( 5,'GAS$      ', 5, 0, 0, 0., 0., STEMAX,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
*
      CALL GSTMED( 6,'BRASS$    ', 6, 0, 0, 0., 0., STEMAX,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
*
      CALL GSTPAR(5,'STRA',1.)
      CALL GSTPAR(5,'DCUTE',5E-5)
*
*             Defines geometry of the set-up
*
*             Basic parameters
      UTHICK = 0.4
      CSHIEL = 0.1
      ESHIEL = 0.1
      WGS1   = 0.55
      WGS2   = 0.35
      BTHICK = 0.05
      AIRGP1 = 0.9
      CTHICK = 2.0
      AIRGP2 = 1.9
      EPCHAM = 1.6
      AIRGP3 = 1.7
      HMODL1 = (2.*CSHIEL+UTHICK+AIRGP1)*.5
      HSECT1 = HMODL1*64.
      HMODL2 = (2.*CSHIEL+UTHICK+AIRGP2)*.5
      HSECT2 = HMODL2*35.
      HMODL3 = (          CTHICK+AIRGP3)*.5
      HSECT3 = HMODL3*13.
      HLCALO = HSECT1+HSECT2+HSECT3+AIRGP1
*
*             Define the overall calorimeter
      PAR(1)=25.
      PAR(2)=25.
      PAR(3)=HLCALO
      CALL GSVOLU('CALO', 'BOX ', 1, PAR, 3, IVOLU)
*
*             Now define the three sections
      PAR(3)=HSECT1
      CALL GSVOLU('CAL1', 'BOX ', 1, PAR, 3, IVOLU)
      PAR(3)=HSECT2
      CALL GSVOLU('CAL2', 'BOX ', 1, PAR, 3, IVOLU)
      PAR(3)=HSECT3
      CALL GSVOLU('CAL3', 'BOX ', 1, PAR, 3, IVOLU)
*
*             and position them
      CALL GSPOS('CAL1',1,'CALO',0.0,0.0,
     * 0.5*AIRGP1-(HSECT2+HSECT3) ,0,'ONLY')
      CALL GSPOS('CAL2',1,'CALO',0.0,0.0,
     * 0.5*AIRGP1+HSECT1-HSECT3  ,0,'ONLY')
      CALL GSPOS('CAL3',1,'CALO',0.0,0.0,
     * 0.5*AIRGP1+HSECT1+HSECT2  ,0,'ONLY')
*
*             Now divide each section in modules
      CALL GSDVN('MOD1','CAL1',64,3)
      CALL GSDVN('MOD2','CAL2',35,3)
      CALL GSDVN('MOD3','CAL3',13,3)
*
*             Define copper shielding for type 1 and 2 modules
      PAR(3)=CSHIEL*.5
      CALL GSVOLU('SHIL', 'BOX ', 2, PAR, 3, IVOLU)
*
*             Define Uranium plate for type 1 and 2 modules
      PAR(3)=UTHICK*.5
      CALL GSVOLU('URPL', 'BOX ', 3, PAR, 3, IVOLU)
*
*             Define brass chamber for type 1 modules
      PAR(3)=(WGS1+BTHICK)*.5
      CALL GSVOLU('CHA1', 'BOX ', 6, PAR, 3, IVOLU)
      CALL GSDVN('TUB1','CHA1', 40, 2)
*
*             Define gas for chamber of type 1 modules
      PAR(1)=24.
      PAR(2)=1.2*.5
      PAR(3)=WGS1*.5
      CALL GSVOLU('GAS1', 'BOX ', 5, PAR, 3, IVOLU)
      CALL GSPOS('GAS1', 1, 'TUB1', 0.0, 0.0, 0.0, 0, 'ONLY')
*
*             Define epoxy wrapping for type 1 chambers
      PAR(1)=25.
      PAR(2)=25.
      PAR(3)=ESHIEL*.5
      CALL GSVOLU('EPO1', 'BOX ', 4, PAR, 3, IVOLU)
*
*             Define copper plate for type 3 modules
      PAR(3)=CTHICK*.5
      CALL GSVOLU('COPL', 'BOX ', 2, PAR, 3, IVOLU)
*
*             Define epoxy chamber for type 2 and 3 modules
      PAR(2)=47.*.5
      PAR(3)=EPCHAM*.5
      CALL GSVOLU('CHA2', 'BOX ', 4, PAR, 3, IVOLU)
      CALL GSDVN('TUB2', 'CHA2', 72, 2)
*
*             Define gas for chamber of type 2 and 3 modules
      PAR(1)=24.
      PAR(2)=(.65-.1)*.5
      PAR(3)=WGS2*.5
      CALL GSVOLU('GAS2', 'BOX ', 5, PAR, 3, IVOLU)
      CALL GSPOS('GAS2', 1, 'TUB2', 0.0, 0.0, 0.0, 0, 'ONLY')
*
*             Now position front chamber
      Z=-HLCALO+0.5*AIRGP1-(ESHIEL+WGS1+BTHICK)*0.5
      CALL GSPOS('EPO1',1,'CALO',0.0,0.0,Z,0,'ONLY')
      Z=-HLCALO+0.5*AIRGP1
      CALL GSPOS('CHA1',1,'CALO',0.0,0.0,Z,0,'ONLY')
      Z=-HLCALO+0.5*AIRGP1+(ESHIEL+WGS1+BTHICK)*0.5
      CALL GSPOS('EPO1',2,'CALO',0.0,0.0,Z,0,'ONLY')
*
*             Now assemble type 1 modules
      Z=CSHIEL*.5-HMODL1
      CALL GSPOS('SHIL',1,'MOD1',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL+UTHICK*.5-HMODL1
      CALL GSPOS('URPL',1,'MOD1',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL*1.5+UTHICK-HMODL1
      CALL GSPOS('SHIL',2,'MOD1',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL*2.+UTHICK+AIRGP1*.5-(WGS1+BTHICK+ESHIEL)*.5-HMODL1
      CALL GSPOS('EPO1',1,'MOD1',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL*2.+UTHICK+AIRGP1*.5-HMODL1
      CALL GSPOS('CHA1',1,'MOD1',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL*2.+UTHICK+AIRGP1*.5+(WGS1+BTHICK+ESHIEL)*.5-HMODL1
      CALL GSPOS('EPO1',2,'MOD1',0.0,0.0,Z,0,'ONLY')
*
*             Now assemble type 2 modules
      Z=CSHIEL*.5-HMODL2
      CALL GSPOS('SHIL',1,'MOD2',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL+UTHICK*.5-HMODL2
      CALL GSPOS('URPL',1,'MOD2',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL*1.5+UTHICK-HMODL2
      CALL GSPOS('SHIL',2,'MOD2',0.0,0.0,Z,0,'ONLY')
      Z=CSHIEL*2.+UTHICK+AIRGP2*.5-HMODL2
      CALL GSPOS('CHA2',1,'MOD2',0.0,0.0,Z,0,'ONLY')
*
*             Now assemble type 3 modules
      Z=CTHICK*.5-HMODL3
      CALL GSPOS('COPL',1,'MOD3',0.0,0.0,Z,0,'ONLY')
      Z=CTHICK+AIRGP3*.5-HMODL3
      CALL GSPOS('CHA2',1,'MOD3',0.0,0.0,Z,0,'ONLY')
*
*             Define geometry optimization
      CALL GSORD('CALO',3)
      CALL GSORD('MOD1',3)
      CALL GSORD('MOD2',3)
      CALL GSORD('MOD3',3)
*
*             Close geometry banks. Mandatory system routine.
*
      CALL GGCLOS
*
      END


      SUBROUTINE UGINIT
*
************************************************************************
*                                                                      *
*          To initialise GEANT3 program and read data cards            *
*                                                                      *
************************************************************************
*
#include "gcbank.inc"  
#include "gcflag.inc"  
#include "gcunit.inc"  
*KEEP,URGEOM.
*
      COMMON/URGEOM/UTHICK,CSHIEL,ESHIEL,BTHICK,AIRGP1,CTHICK,AIRGP2,
     *              AIRGP3,EPCHAM,HMODL1,HSECT1,HMODL2,HSECT2,HMODL3,
     *              HSECT3,HLCALO,WGS1,WGS2,IGAST
C
*KEEP,URMIPS.
      COMMON/URMIPS/XMIP,NTXMIP,TXMIP,TXMIP2,CALIB1,CALIB2,XKPL(60,8),
     *              XPHMIP(20,34),IDSTW
C
*KEND.
*
*     -----------------------------------------------------------------
*
*             Open user files
*
      CALL UFILES
*
*             Initialize GEXAM4 global variables
*
      IGAST = 1
      IDSTW = 0
      NTXMIP= 0
      TXMIP = 0.
      TXMIP2= 0.
      CALL VZERO(XKPL,300)
*
*             Initialize GEANT
      CALL GINIT
*
*             Prints version number
      WRITE(LOUT,1001)
*
*             Define a data card 'GAST' to change gas type
*
      CALL FFKEY('GAST',IGAST,1,'INTEGER')
*
*             Define a data card to fill and write the PH matrix
*
      CALL FFKEY('DSTW',IDSTW,1,'INTEGER')
*
*             Read data cards with FFREAD
*
      write(lout,1000)
 1000 format(/,' ========> Reading ffread data cards : type <======='
     +,/,'read 4'
     +,/,'your own data cards if any'
     +,/,'stop',/,'      Now waiting for input',/)

      CALL GFFGO
*
*             Initialize GEANT/ZBOOK data structures
      CALL GZINIT
*
*
*             Geometry and materials description.
      CALL UGEOM
*
*             Particle table definition and energy loss initialization.
      CALL GPART
      CALL GPHYSI
*
*             Compute calibration factors for the 2 types of chambers.
*             Calibration is obtained with minimum ionizing muons
*             of 500 MeV in material 5 (Argon/Isobuthane or Argon/CO2)
*
      JMA    = LQ(JMATE-5)
      JMULOS = LQ(JMA-2)
      DEDX   = Q(JMULOS+45)
      CALIB1 = DEDX*WGS1
      CALIB2 = DEDX*WGS2
*
      CALL UHINIT
*
 1001 FORMAT(/,' **** GEXAM4 VERSION 1.17.00 ( 10 May 1988 ) ',/)
*
      END


      SUBROUTINE UGLAST
*
************************************************************************
*                                                                      *
*       Termination routine to print histograms and statistics         *
*                                                                      *
************************************************************************
*
*KEEP,URMIPS.
      COMMON/URMIPS/XMIP,NTXMIP,TXMIP,TXMIP2,CALIB1,CALIB2,XKPL(60,8),
     *              XPHMIP(20,34),IDSTW
C
*KEND.
*
*     -----------------------------------------------------------------
*
      CALL GLAST
*
      CALL PXMIPS
*
      CALL HNOENT(10,NOENT)
C
      DO 20 I=1,5
         IF(NOENT.GT.0)THEN
            DO 10 J=1,60
              XKPL(J,I)=XKPL(J,I)/FLOAT(NOENT)
  10        CONTINUE
         ENDIF
         CALL HPAK(100+I,XKPL(1,I))
  20  CONTINUE
C
C             Save histograms
C
      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('HBOOK')
*
*             Print HBOOK histograms
      CALL HPRINT(0)
*
      END


      SUBROUTINE UHINIT
*
************************************************************************
*                                                                      *
*             To book the user's histograms                            *
*                                                                      *
************************************************************************
*
#include "gckine.inc"  

      PARAMETER (NKEYS=6)
      CHARACTER*8 CHNAME(NKEYS)
      DATA CHNAME/'X','Y','Z','GEKIN','TOFG','FLAG'/
*
*     ------------------------------------------------------------------
*
      XMIN=0.
      XMAX=100.*PKINE(1)
      CALL HBOOKN(1,'NTUPLE',NKEYS,'//HBOOK',9950,CHNAME)
      CALL HBOOK1(5,'PULSE HEIGHT DISTRIBUTION IN MIPS$',100,0.,100.,0.)
      CALL HBOOK1(6,'PULSE HEIGHT DISTRIBUTION IN MIPS$',100,0.,2.,0.)
      CALL HBOOK1(10,'NUMBER OF MIPS PER EVENT$',100,XMIN,XMAX,0.)
      CALL HBOOK1(101,'NUMBER OF PHOTONS PER PLANE$',60,1.,60.,0.)
      CALL HCOPY(101,102,'NUMBER OF ELECTRONS+POSITRONS PER PLANE$')
      CALL HCOPY(101,103,'NUMBER OF NEUTRAL PER PLANE$')
      CALL HCOPY(101,104,'NUMBER OF CHARGED HADRONS PER PLANE$')
      CALL HCOPY(101,105,'NUMBER OF MUONS PER PLANE$')
      CALL HBOOK1(500,'MUON AVERAGE ENERGY LOSS$',100,1.E-6,5.E-6,0.)
*
      END


      SUBROUTINE GPHXSD
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Debugs  PHXS bank containing x-section constants              *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gcbank.inc"  
*KEEP,GCJLOC.
      COMMON/GCJLOC/NJLOC(2),JTM,JMA,JLOSS,JPROB,JMIXT,JPHOT,JANNI
     +                  ,JCOMP,JBREM,JPAIR,JDRAY,JPFIS,JMUNU,JRAYL
     +                  ,JMULOF,JCOEF,JRANG
C
      INTEGER       NJLOC   ,JTM,JMA,JLOSS,JPROB,JMIXT,JPHOT,JANNI
     +                  ,JCOMP,JBREM,JPAIR,JDRAY,JPFIS,JMUNU,JRAYL
     +                  ,JMULOF,JCOEF,JRANG
C
      COMMON/GCJLCK/NJLCK(2),JTCKOV,JABSCO,JEFFIC,JINDEX,JCURIN
     +                      ,JPOLAR,JTSTRA,JTSTCO,JTSTEN,JTASHO
C
      EQUIVALENCE (JLASTV,JTSTEN)
C
      INTEGER       NJLCK,JTCKOV,JABSCO,JEFFIC,JINDEX,JCURIN
     +                   ,JPOLAR,JLASTV,JTSTRA,JTSTCO,JTSTEN
     +                   ,JTASHO
C
*KEEP,GCPMXZ.
      INTEGER MAXELZ
      PARAMETER (MAXELZ=100)
C
*KEEP,GCPHXS.
      INTEGER MAXPOW,MAXINT
      PARAMETER (MAXPOW=4)
      PARAMETER (MAXINT=13)
      CHARACTER*6 CRNGUP
      COMMON /GCPXRN/ CRNGUP(MAXINT,MAXELZ)
      REAL COFS,GPOMIN
      COMMON /GCPXCF/ COFS(MAXPOW,MAXINT,MAXELZ),GPOMIN(MAXELZ)
C
*KEEP,GC10EV.
      REAL G10EV,TENEV
      PARAMETER (G10EV=1.0E-8)
      PARAMETER (TENEV=1.E-2)
C
*KEND.
      CHARACTER*20 CHTMED, CHMATE
*
* Loop over material
      DO 10 ITM=1,IQ(JTMED-2)
         JTM=LQ(JTMED-ITM)
         IF(JTM.GT.0) THEN
            IMA=Q(JTM+6)
            JMA=LQ(JMATE-IMA)
            CALL UHTOC(IQ(JTM+1),4,CHTMED,20)
            CALL UHTOC(IQ(JMA+1),4,CHMATE,20)
            WRITE(6,10000) CHTMED(:LNBLNK(CHTMED)),
     +      CHMATE(:LNBLNK(CHMATE))
            JPHOT = LQ(JMA-6)
            JPHXS=LQ(JPHOT-1)
            NZ=Q(JPHXS+1)
            NIT=Q(JPHXS+3*NZ+2)
            WRITE(6,10100) NZ
            WRITE(6,10200) (Q(JPHXS+1+IZ),Q(JPHXS+1+2*NZ+IZ),IZ=1,NZ)
            WRITE(6,10300)(IT+1,(Q(JPHXS+3*NZ+2+5*IT+IP),IP=1,5),IT=0,
     +      NIT- 1)
         ENDIF
   10 CONTINUE
10000 FORMAT(' Tracking medium : ',A,' Material : ',A)
10100 FORMAT(1X,I3,' elements in this material with ',I3,
     +       ' energy intervals')
10200 FORMAT('    Z    = ',F5.2,'    Weight   = ',F5.2)
10300 FORMAT((' Int. N. ',I2,' E= ',E11.4,' K= ',4(E11.4,',')))
      END
