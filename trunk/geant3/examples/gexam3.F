      SUBROUTINE UFILES
*
*            To open FFREAD and HBOOK files
*
      CHARACTER*(*) FILNAM, FSTAT
      PARAMETER (FILNAM='gexam3.dat')
*
      PARAMETER (FSTAT='OLD')
*
      OPEN(UNIT=4,FILE=FILNAM,STATUS=FSTAT,FORM='FORMATTED')
*
*             Open a HBOOK4 direct access file
*
      CALL HROPEN(34,'HBOOK','gexam3.hist','N',1024,ISTAT)
      END


      SUBROUTINE UGINIT
C
#include "gclist.inc"
*KEEP,GCRZ.
      COMMON/GCRZ1/NRECRZ,NRGET,NRSAVE,LRGET(20),LRSAVE(20)
      INTEGER      NRECRZ,NRGET,NRSAVE,LRGET    ,LRSAVE
      COMMON/GCRZ2/RZTAGS
      CHARACTER*8 RZTAGS(4)
C
*KEND.
C
      COMMON/DLSFLD/ISWFLD,FLDVAL
C
      DIMENSION AMYL(3),ZMYL(3),WMYL(3)
      DIMENSION APYR(4),ZPYR(4),WPYR(4)
      DIMENSION AFRE(3),ZFRE(3),WFRE(3)
      DIMENSION ACO2(2),ZCO2(2),WCO2(2)
      DIMENSION APOL(2),ZPOL(2),WPOL(2)
      DIMENSION ASCI(2),ZSCI(2),WSCI(2)
      DIMENSION AFDR(5),ZFDR(5),WFDR(5)
      DIMENSION ARDR(5),ZRDR(5),WRDR(5)
      CHARACTER*4 CHGET(20),CHSAVE(20)
C
      DATA AMYL,ZMYL,WMYL/12.01,1.01,16.00,6.,1.,8.,5.,4.,2./
      DATA APYR,ZPYR,WPYR/28.09,10.82,23.00,16.00,14.,5.,11.,8.,
     + .385,.039,.038,.538/
      DATA AFRE,ZFRE,WFRE/12.01,35.45,19.00,6.,17.,9.,1.,1.,3./
      DATA ACO2,ZCO2,WCO2/12.01,16.00,6.,8.,1.,2./
      DATA APOL,ZPOL,WPOL/12.01,1.01,6.,1.,1.,2./
      DATA ASCI,ZSCI,WSCI/12.01,1.01,6.,1.,1.,1./
      DATA AFDR,ZFDR,WFDR/63.54,183.85,39.95,12.01,1.01,
     + 29.,74.,18.,6.,1.,.457,.009,.313,.177,.045/
      DATA ARDR,ZRDR,WRDR/63.54,183.85,39.95,12.01,1.01,
     + 29.,74.,18.,6.,1.,.528,.002,.275,.156,.039/
*
*             Open user files
*
      CALL UFILES
C
      CALL GINIT
C
C                 Define here user data cards
C
      ISWFLD=1
      FLDVAL=3.9
      CALL FFKEY('FIEL',ISWFLD,2,'MIXED')
C
      write(6,1000)
 1000 format(/,' ========> Reading ffread data cards : type <======='
     +,/,'read 4'
     +,/,'your own data cards if any'
     +,/,'stop',/,'      Now waiting for input',/)

      CALL GFFGO
C
      CALL GZINIT
C
      CALL GDINIT
C
      IF(NRGET.NE.0)THEN
         IF(IUCOMP(4HINIT,LRGET,NRGET).NE.0)THEN
            CALL GRFILE(3,'TEST',' ')
            CALL GRGET('INIT',0,1,0)
            GO TO 10
         ENDIF
      ENDIF
C
      IF(NGET .GT.0) CALL GOPEN(1,'I',0,IER)
      IF(NSAVE.GT.0) CALL GOPEN(2,'O',0,IER)
C
      DO 5 I=1,20
         CALL UHTOC(LGET(I),4,CHGET(I),4)
         CALL UHTOC(LSAVE(I),4,CHSAVE(I),4)
   5  CONTINUE
      CALL GGET(1,CHGET,-NGET,IDENT,IER)
      IF(IDENT.EQ.0) GO TO 10
C
      CALL GPART
C
      CALL GSMATE(15,'AIR$',14.61,7.3,0.001205,30423.24,6750.,0,0)
      CALL GSMATE(16,'VACUUM$',1.E-16,1.E-16,1.E-16,1.E+16,1.E+16,0,0)
      CALL GSMATE(17,'CALCIUM$',40.08,20.,1.55,10.4,23.2,0,0)
      CALL GSMIXT(18,'MYLAR$',AMYL,ZMYL,1.39,-3,WMYL)
      CALL GSMIXT(19,'PYREX$',APYR,ZPYR,2.23,4,WPYR)
      CALL GSMIXT(20,'FREON13$',AFRE,ZFRE,4.26E-03,-3,WFRE)
      CALL GSMIXT(21,'CO2$',ACO2,ZCO2,1.79E-03,-2,WCO2)
      CALL GSMIXT(22,'POLYETHYLENE$',APOL,ZPOL,.935,-2,WPOL)
      CALL GSMIXT(23,'SCINT$',ASCI,ZSCI,1.032,-2,WSCI)
      CALL GSMIXT(24,'FRONTDRIFT$',AFDR,ZFDR,2.85E-03,5,WFDR)
      CALL GSMIXT(25,'REARDRIFT$',ARDR,ZRDR,3.24E-03,5,WRDR)
C
      CALL UGEOM
C
      CALL UDET
C
  10  CALL GPHYSI
C
      IF(NHSTA.GT.0) CALL GBHSTA
      CALL USETPL
C
C
      RETURN
      END


      SUBROUTINE UGLAST
C
#include "gclist.inc"
C
      CALL GLAST
C
      IF(NGET.NE.0.OR.NSAVE.NE.0) CALL GCLOSE(0,IER)
C
      CALL HISTDO
C
      RETURN
      END


      SUBROUTINE UGEOM
C
C
#include "gclist.inc"
#include "gconst.inc"
C
      COMMON/DLSFLD/ISWFLD,FLDVAL
C
      DIMENSION   CAVPAR(3),TGTPAR(3),TBIPAR(4),TBOPAR(4),ARMPAR(4)
     + ,FMIPAR(3),FCIPAR(5),FCOPAR(5),FHOPAR(3),FDIPAR(3),FDOPAR(3)
     + ,FLDPAR(3),RD1IPA(3),RD1OPA(3),RD2IPA(3),RD2OPA(3),RD3IPA(3)
     + ,RD3OPA(3),RHOPAR(3)
C
C      VOLUME SIZE PARAMETERS
C
      DATA CAVPAR/450.,100.,360./
      DATA TGTPAR/.47625,1.27,.0232/
      DATA TBIPAR/60.91172,5.07172,10.14,27.92/
      DATA TBOPAR/60.96,5.08,10.16,27.94/
      DATA ARMPAR/18.,210.,50.,190./
      DATA FMIPAR/36.83,12.7,.06/
      DATA FCIPAR/13.29941,50.69885,7.565299,13.90419,34.23/
      DATA FCOPAR/13.335,50.80,7.62,13.97,34.29/
      DATA FHOPAR/40.005,13.49375,.15875/
      DATA FDIPAR/41.91,13.97,4.2/
      DATA FDOPAR/41.917,13.977,4.207/
      DATA FLDPAR/55.88,19.05,25.40/
      DATA RD1IPA/70.000,25.500,7.800/
      DATA RD1OPA/70.007,25.507,7.807/
      DATA RD2IPA/85.000,29.500,7.800/
      DATA RD2OPA/85.007,29.507,7.807/
      DATA RD3IPA/100.000,33.500,7.800/
      DATA RD3OPA/100.007,33.507,7.807/
      DATA RHOPAR/150.5,47.,.3175/
C
C      LOCATIONS:
C          TGT VERSUS CAVE CENTER
C          FIELD AND DETECTORS VERSUS ARM CENTER
C          ANGLE SETTINGS OF LEFT AND RIGHT ARMS (DEG., ANGR IS NEG.)
C          LEFT AND RIGHT ARM VERSUS CAVE CENTER
C
      DATA ZTG/-100./
      DATA DFCO,DFHO,DFDO,DFLD,DRD1,DRD2,DRD3,DRHO/-154.77,-120.1625
     + ,-115.4,-84.92,-27.77,3.98,35.73,142.41/
      DATA ANGL,ANGR/45.,-45./
      DATA DARM/217./
      CANGL=COS(ANGL*DEGRAD)
      SANGL=SIN(ANGL*DEGRAD)
      CANGR=COS(ANGR*DEGRAD)
      SANGR=SIN(ANGR*DEGRAD)
      XLARM=DARM*SANGL
      ZLARM=DARM*CANGL+ZTG
      XRARM=DARM*SANGR
      ZRARM=DARM*CANGR+ZTG
C
      CALL GSTMED( 1,'AIR       $',15,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED( 2,'TARGET    $',17,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED( 3,'VACUUM    $',16,0,0,0., 0.,.0,.0,.001,.5,0,0)
      CALL GSTMED( 4,'MYLAR     $',18,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED( 5,'PYREX     $',19,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED( 6,'FREON13   $',20,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED( 7,'POLYETHYL.$',22,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED( 8,'SCINTILL. $',23,0,0,0.,10.,.2,.1,.01,.5,0,0)
      CALL GSTMED( 9,'FRONTDRIFT$',24,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED(10,'AIR+FIELD $',15,0,ISWFLD,6.,10.,.2,.1,.01,.5,0,0)
      CALL GSTMED(11,'REAR DRIFT$',25,0,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED(12,'CO2       $',21,0,0,0.,10.,.2,.1,.01,.5,0,0)
      CALL GSTMED(13,'FR.DR.S.P.$',24,1,0,0.,10.,.2,.1,.001,.5,0,0)
      CALL GSTMED(14,'RE.DR.S.P.$',25,1,0,0.,10.,.2,.1,.001,.5,0,0)
C
C      CAVE REFERENCE FRAME:
C          OZ ALONG BEAM, OX HORIZONTAL, OY VERTICAL
C
      CALL GSVOLU('CAVE','BOX ',1,CAVPAR,3,ICAVE)
C
C      TARGET BOX SHIFTED UPSTREAM 100CM IN CAVE.
C           REFERENCE FRAME:
C               OZ ALONG BEAM, OX HORIZONTAL, OY VERTICAL.
C           CA DENSITY = 1.55 G/CM**3.
C           5 SEGMENTS.
C
      CALL GSVOLU('TGT ','BOX ',2,TGTPAR,3,ITGT )
      CALL GSVOLU('TBIN','TRD1',3,TBIPAR,4,ITBIN)
      CALL GSVOLU('TBOU','TRD1',4,TBOPAR,4,ITBOU)
C
C      LEFT AND RIGHT ARM REFERENCE FRAMES:
C           OZ ALONG CENTRAL RAY, OX HORIZONTAL, OY VERTICAL.
C
      CALL GSVOLU('ARM ','TRD1',1,ARMPAR,4,IARM)
      CALL GSVOLU('FMIR','BOX ',5,FMIPAR,3,IFMIR)
      CALL GSVOLU('FCIN','TRD2',6,FCIPAR,5,IFCIN)
      CALL GSVOLU('FCOU','TRD2',7,FCOPAR,5,IFCOU)
      CALL GSVOLU('FHOD','BOX ',8,FHOPAR,3,IFHOD)
      CALL GSVOLU('FDIN','BOX ',9,FDIPAR,3,IFDIN)
      CALL GSVOLU('FDOU','BOX ',4,FDOPAR,3,IFDOU)
      CALL GSVOLU('FLD ','BOX ',10,FLDPAR,3,IFLD )
      CALL GSVOLU('RD1I','BOX ',11,RD1IPA,3,IRD1I)
      CALL GSVOLU('RD1O','BOX ', 4,RD1OPA,3,IRD1O)
      CALL GSVOLU('RD2I','BOX ',11,RD2IPA,3,IRD2I)
      CALL GSVOLU('RD2O','BOX ', 4,RD2OPA,3,IRD2O)
      CALL GSVOLU('RD3I','BOX ',11,RD3IPA,3,IRD3I)
      CALL GSVOLU('RD3O','BOX ', 4,RD3OPA,3,IRD3O)
      CALL GSVOLU('RHOD','BOX ', 8,RHOPAR,3,IRHOD)
C
C         DRIFT CHAMBER SENSE PLANES
C
      FDIPAR(3)=.0025
      CALL GSVOLU('FSP ','BOX ',13,FDIPAR,3,IFSP )
      RD1IPA(3)=.0025
      CALL GSVOLU('RSP1','BOX ',14,RD1IPA,3,IRSP1)
      RD2IPA(3)=.0025
      CALL GSVOLU('RSP2','BOX ',14,RD2IPA,3,IRSP2)
      RD3IPA(3)=.0025
      CALL GSVOLU('RSP3','BOX ',14,RD3IPA,3,IRSP3)
C
C
      THLX=90.+ANGL
      THRX=90.+ANGR
      PHX=0.
      THY=90.
      PHY=90.
      THLZ=ANGL
      THRZ=-ANGR
      PHLZ=0.
      PHRZ=180.
      CALL GSROTM(1,THLX,PHX,THY,PHY,THLZ,PHLZ)
      CALL GSROTM(2,THRX,PHX,THY,PHY,THRZ,PHRZ)
C
C
      CALL GSPOS('TGT ',1,'TBIN', 0., 0.,-5.08,0,'ONLY')
      CALL GSPOS('TGT ',2,'TBIN', 0., 0.,-2.54,0,'ONLY')
      CALL GSPOS('TGT ',3,'TBIN', 0., 0., 0.  ,0,'ONLY')
      CALL GSPOS('TGT ',4,'TBIN', 0., 0., 2.54,0,'ONLY')
      CALL GSPOS('TGT ',5,'TBIN', 0., 0., 5.08,0,'ONLY')
      CALL GSPOS('TBIN',1,'TBOU', 0., 0.,   0.,0,'ONLY')
      CALL GSPOS('TBOU',1,'CAVE', 0., 0.,  ZTG,0,'ONLY')
C
      CALL GSPOS('ARM ',1,'CAVE',XLARM,0.,ZLARM,1,'ONLY')
      CALL GSPOS('ARM ',2,'CAVE',XRARM,0.,ZRARM,2,'ONLY')
C
      CALL GSPOS('FMIR',1,'FCIN',0.,0.,26.0350,0,'ONLY')
      CALL GSPOS('FCIN',1,'FCOU',0.,0., 0.    ,0,'ONLY')
      CALL GSPOS('FCOU',1,'ARM ',0.,0., DFCO  ,0,'ONLY')
C
      CALL GSPOS('FHOD',1,'ARM ',0.,0., DFHO  ,0,'ONLY')
C
      CALL GSPOS('FDIN',1,'FDOU',0.,0., 0.    ,0,'ONLY')
      CALL GSPOS('FDOU',1,'ARM ',0.,0., DFDO  ,0,'ONLY')
C
      CALL GSPOS('FLD ',1,'ARM ',0.,0., DFLD  ,0,'ONLY')
C
      CALL GSPOS('RD1I',1,'RD1O',0.,0., 0.    ,0,'ONLY')
      CALL GSPOS('RD1O',1,'ARM ',0.,0., DRD1  ,0,'ONLY')
C
      CALL GSPOS('RD2I',1,'RD2O',0.,0., 0.    ,0,'ONLY')
      CALL GSPOS('RD2O',1,'ARM ',0.,0., DRD2  ,0,'ONLY')
C
      CALL GSPOS('RD3I',1,'RD3O',0.,0., 0.    ,0,'ONLY')
      CALL GSPOS('RD3O',1,'ARM ',0.,0., DRD3  ,0,'ONLY')
C
      CALL GSPOS('RHOD',1,'ARM ',0.,0., DRHO  ,0,'ONLY')
C
      CALL GSPOS('FSP ',1,'FDIN',0.,0.,-2.9975,0,'ONLY')
      CALL GSPOS('FSP ',2,'FDIN',0.,0.,-1.7975,0,'ONLY')
      CALL GSPOS('FSP ',3,'FDIN',0.,0., -.5975,0,'ONLY')
      CALL GSPOS('FSP ',4,'FDIN',0.,0.,  .6025,0,'ONLY')
      CALL GSPOS('FSP ',5,'FDIN',0.,0., 1.8025,0,'ONLY')
      CALL GSPOS('FSP ',6,'FDIN',0.,0., 3.0025,0,'ONLY')
C
      CALL GSPOS('RSP1',1,'RD1I',0.,0.,-5.9975,0,'ONLY')
      CALL GSPOS('RSP1',2,'RD1I',0.,0.,-3.5975,0,'ONLY')
      CALL GSPOS('RSP1',3,'RD1I',0.,0.,-1.1975,0,'ONLY')
      CALL GSPOS('RSP1',4,'RD1I',0.,0., 1.2025,0,'ONLY')
      CALL GSPOS('RSP1',5,'RD1I',0.,0., 3.6025,0,'ONLY')
      CALL GSPOS('RSP1',6,'RD1I',0.,0., 6.0025,0,'ONLY')
C
      CALL GSPOS('RSP2',1,'RD2I',0.,0.,-5.9975,0,'ONLY')
      CALL GSPOS('RSP2',2,'RD2I',0.,0.,-3.5975,0,'ONLY')
      CALL GSPOS('RSP2',3,'RD2I',0.,0.,-1.1975,0,'ONLY')
      CALL GSPOS('RSP2',4,'RD2I',0.,0., 1.2025,0,'ONLY')
      CALL GSPOS('RSP2',5,'RD2I',0.,0., 3.6025,0,'ONLY')
      CALL GSPOS('RSP2',6,'RD2I',0.,0., 6.0025,0,'ONLY')
C
      CALL GSPOS('RSP3',1,'RD3I',0.,0.,-5.9975,0,'ONLY')
      CALL GSPOS('RSP3',2,'RD3I',0.,0.,-3.5975,0,'ONLY')
      CALL GSPOS('RSP3',3,'RD3I',0.,0.,-1.1975,0,'ONLY')
      CALL GSPOS('RSP3',4,'RD3I',0.,0., 1.2025,0,'ONLY')
      CALL GSPOS('RSP3',5,'RD3I',0.,0., 3.6025,0,'ONLY')
      CALL GSPOS('RSP3',6,'RD3I',0.,0., 6.0025,0,'ONLY')
C
C      CALL GSORD('ARM ',3)
      CALL GSORD('FDIN',3)
      CALL GSORD('RD1I',3)
      CALL GSORD('RD2I',3)
      CALL GSORD('RD3I',3)
C
C
      IF(IUCOMP(4HVOLU,LPRIN,NPRIN).NE.0)CALL GPVOLU(0)
      IF(IUCOMP(4HROTM,LPRIN,NPRIN).NE.0)CALL GPROTM(0)
      IF(IUCOMP(4HTMED,LPRIN,NPRIN).NE.0)CALL GPTMED(0)
      IF(IUCOMP(4HMATE,LPRIN,NPRIN).NE.0)CALL GPMATE(0)
      IF(IUCOMP(4HPART,LPRIN,NPRIN).NE.0)CALL GPPART(0)
C
C
      CALL GGCLOS
C
      RETURN
      END


      SUBROUTINE UDET
C
#include "gclist.inc"
      CHARACTER*4 NAMESH(9),NAFD(2),NARD1(2),NARD2(2),NARD3(2)
      DIMENSION NBITSH(9),NBITSV(2),ORIG(9),FACT(9)
C
      DATA NAFD,NARD1,NARD2,NARD3/'ARM ','FSP ','ARM ','RSP1'
     + ,'ARM ','RSP2','ARM ','RSP3'/
      DATA NBITSV/2,6/
      DATA NAMESH/'X   ','Y   ','Z   ','U   ','V   ','W   '
     + ,'UP  ','VP  ','WP  '/
      DATA NBITSH/6*17,3*12/
      DATA ORIG/3*300,3*100.,3*1./
      DATA FACT/3*100.,3*500.,3*1000./
C
C
      IF(IUCOMP(4HDRFT,LSETS,NSETS).EQ.0)GO TO 99
      CALL GSDET('DRFT','FSP ',2,NAFD ,NBITSV,1,100,100,IDRFT,IFD )
      CALL GSDET('DRFT','RSP1',2,NARD1,NBITSV,1,100,100,IDRFT,IRD1)
      CALL GSDET('DRFT','RSP2',2,NARD2,NBITSV,1,100,100,IDRFT,IRD2)
      CALL GSDET('DRFT','RSP3',2,NARD3,NBITSV,1,100,100,IDRFT,IRD3)
C
C
      CALL GSDETH('DRFT','FSP ',9,NAMESH,NBITSH,ORIG,FACT)
      CALL GSDETH('DRFT','RSP1',9,NAMESH,NBITSH,ORIG,FACT)
      CALL GSDETH('DRFT','RSP2',9,NAMESH,NBITSH,ORIG,FACT)
      CALL GSDETH('DRFT','RSP3',9,NAMESH,NBITSH,ORIG,FACT)
C
C
      IF(IUCOMP(4HSETS,LPRIN,NPRIN).NE.0)CALL GPSETS(0,0)
C
  99  RETURN
      END


      SUBROUTINE USETPL
C
#include "gclist.inc"
C
      CALL HBOOK1(100,' XHITS IN LFDR$',100,-50.,50.,0.)
      CALL HBOOK1(101,' YHITS IN LFDR$', 40,-20.,20.,0.)
      CALL HBOOK1(102,' ZHITS IN LFDR$', 20,-10.,10.,0.)
      CALL HBOOK1(110,' XHITS IN RFDR$',100,-50.,50.,0.)
      CALL HBOOK1(111,' YHITS IN RFDR$', 40,-20.,20.,0.)
      CALL HBOOK1(112,' ZHITS IN RFDR$', 20,-10.,10.,0.)
C
      RETURN
      END


      SUBROUTINE GUKINE
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       Read or Generates Kinematics for primary tracks          *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gclist.inc"
#include "gckine.inc"
#include "gconst.inc"
      CHARACTER*4 CHGET(20)
      DIMENSION XVERT(3)
      DIMENSION P1(3)
      DIMENSION RNDM(3)
      SAVE XVERT
C
      DATA XVERT/0.,0.,-100./
C.
C.    ------------------------------------------------------------------
C.
      IF(NGET.LE.0)GO TO 10
      ITRY=0
  5   ITRY=ITRY+1
      IF(ITRY.GT.2)GO TO 10
      DO 6 I=1,NGET
         CALL UHTOC(LGET(I),4,CHGET(I),4)
   6  CONTINUE
      CALL GGET(1,CHGET,NGET,IDENT,IER)
      IF(IDENT.LE.0)GO TO 5
      IF(IUCOMP(4HKINE,LGET,NGET).NE.0)GO TO 99
C
 10   CALL GSVERT(XVERT,0,0,0,0,NVTX)
C
      IF(IKINE.EQ.0)THEN
           READ(4,100)NTR
  100      FORMAT(I5)
           DO 200 ITR=1,NTR
           READ(4,110)P1,ITYPE
  110      FORMAT(3F8.4,I2)
           CALL GSKINE(P1,ITYPE,NVTX,0,0,NT1)
  200      CONTINUE
      ENDIF
C
      IF(IKINE.EQ.1)THEN
           IF(PKINE(1).EQ.0.)PKINE(1)=25.
           IF(PKINE(2).EQ.0.)PKINE(2)=65.
           IF(PKINE(3).EQ.0.)PKINE(3)=-5.
           IF(PKINE(4).EQ.0.)PKINE(4)= 5.
           IF(PKINE(5).EQ.0.)PKINE(5)=  .2
           IF(PKINE(6).EQ.0.)PKINE(6)=  .8
           IF(PKINE(7).EQ.0.)PKINE(7)= 1.
           IMAX=INT(PKINE(7))
           IF(PKINE(8).EQ.0.)PKINE(8)=14.
           ITYPE=INT(PKINE(8))
           CALL GRNDM(RNDM,3)
           DO 300 I=1,IMAX
           TH=RNDM(1)*(PKINE(2)-PKINE(1))+PKINE(1)
           PH=RNDM(2)*(PKINE(4)-PKINE(3))+PKINE(3)
           PMOM=RNDM(3)*(PKINE(6)-PKINE(5))+PKINE(5)
           P1(1)=PMOM*SIN(TH*DEGRAD)*COS(PH*DEGRAD)
           P1(2)=PMOM*SIN(TH*DEGRAD)*SIN(PH*DEGRAD)
           P1(3)=PMOM*COS(TH*DEGRAD)
           CALL GSKINE(P1,ITYPE,NVTX,0,0,NT1)
  300      CONTINUE
      ENDIF
C
  99  RETURN
      END


      SUBROUTINE GUTREV
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       User routine to control tracking of one event            *
C.    *                                                                *
C.    *       Called by GRUN                                           *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
C.
C.    ------------------------------------------------------------------
C.
#include "gcbank.inc"
#include "gclist.inc"
#include "gcflag.inc"
C
      IF(NGET.EQ.0)GO TO 10
      IF(IUCOMP(4HHITS,LGET,NGET).NE.0)GO TO 99
C
C
  10  CALL GTREVE
C
  99  RETURN
      END


      SUBROUTINE GUFLD(X,F)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       User routine to compute the magnetic field F             *
C.    *       at space point X                                         *
C.    *       Called by GRkuta,GHELIX                                  *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      COMMON/DLSFLD/ISWFLD,FLDVAL
C
#ifndef SINGLEFIELD
      DOUBLE PRECISION X, F
#endif
      DIMENSION X(3),F(3)
C.
C.    ------------------------------------------------------------------
C.
      F(1)=0.
      F(2)=FLDVAL
      F(3)=0.
C
      RETURN
      END


      SUBROUTINE GUSTEP
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       User routine called at the end of each tracking step     *
C.    *       INWVOL is different from 0 when the track has reached    *
C.    *              a volume boundary                                 *
C.    *       ISTOP is different from 0 if the track has stopped       *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gcflag.inc"
#include "gctmed.inc"
#include "gckine.inc"
#include "gctrak.inc"
#include "gcking.inc"
**KEEP,GCSETS.
      COMMON/GCSETS/IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV(20)
C
      INTEGER       IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV
C
*KEND.
      DIMENSION HITS(9)
C.
C.    ------------------------------------------------------------------
C.
c      if (iswit(6).ne.0) call grecord

      IF(IDET.EQ.0)GO TO 20
      IF(CHARGE.EQ.0.)GO TO 20
      IF(INWVOL.NE.1)GO TO 20
      DO 10 I=1,3
         HITS(I) =VECT(I)
   10 CONTINUE
      CALL GMTOD(VECT(1),HITS(4),1)
      CALL GMTOD(VECT(4),HITS(7),2)
      CALL GSAHIT(ISET,IDET,ITRA,NUMBV,HITS,IHIT)
C
   20 IF(ISWIT(6).EQ.0)THEN
C
         IF(NGKINE.GT.0)THEN
            DO 30 I=1,NGKINE
               IFLGK(I)=1
   30       CONTINUE
            CALL GSKING(0)
         ENDIF
      ENDIF
C
      CALL GDEBUG
C
  999 END


      SUBROUTINE GUOUT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       User routine called at the end of each event             *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gclist.inc"
#include "gcflag.inc"
      DIMENSION NUMVS(2),NUMBV(2,50),HITS(9,50),ITRA(50)
      CHARACTER*4 CHSAVE(20)
      SAVE NUMVS
      DATA NUMVS/0,1/
C.
C.    ------------------------------------------------------------------
C.
      IF(IDEBUG.NE.0)THEN
        IF(ISWIT(1).EQ.1)THEN
           CALL GPVERT(0)
           CALL GPKINE(0)
        ENDIF
        IF(ISWIT(2).EQ.1)THEN
           CALL GPJXYZ(0)
        ENDIF
        IF(ISWIT(3).EQ.1)THEN
           CALL GPHITS('*','*')
        ENDIF
      ENDIF
C
      CALL GFHITS('DRFT','FSP ',2,9,50,0,NUMVS,ITRA,NUMBV,HITS,NHITS)
C
      IF(NHITS.LE.0)GO TO 20
      DO 10 IHIT=1,NHITS
      ID=0
      IF(NUMBV(1,IHIT).EQ.1)ID=100
      IF(NUMBV(1,IHIT).EQ.2)ID=110
      X=HITS(4,IHIT)
      Y=HITS(5,IHIT)
      Z=HITS(6,IHIT)
      IF(ID.EQ.0)GO TO 10
      CALL HFILL(ID  ,X,0.,1.)
      CALL HFILL(ID+1,Y,0.,1.)
      CALL HFILL(ID+2,Z,0.,1.)
   10 CONTINUE
C
C                    SAVE EVENT
C
  20  IF(NSAVE.LE.0)GO TO 30
      DO 25 I=1,NSAVE
         CALL UHTOC(LSAVE(I),4,CHSAVE(I),4)
  25  CONTINUE
      CALL GSAVE(2,CHSAVE,NSAVE,1,IER)
C
  30  CALL DLSOUT
C
   99 RETURN
      END


      SUBROUTINE GUDIGI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *       User routine to digitize one event                       *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "gcbank.inc"
#include "gcflag.inc"
#include "gcnum.inc"
C      DIMENSION NUMVS(1),NUMBV(1,30),HITS(7,30),ITRA(30)
C      DIMENSION XL(3),XR(3),PARL(3),PARR(3)
C      DIMENSION X(3),XNEW(3),XP(3),DX(3),XPNEW(3)
C      DIMENSION X1(6),X2(6),XINT(6),PZINT(4)
C      DIMENSION XSAV(2,6,30),YSAV(2,6,30)
C      DATA NUMVS/0/
C      DATA DX/0.,0.,0./
C      DATA NPL/6/
C      DATA S1,S2/0.,0./
C.
C.    ------------------------------------------------------------------
C.
C
C          READING OF TRANSFORMATION AND SHAPE PARAMETERS
C
C      IF(ISWIT(10).EQ.0)GO TO 99
C      IVO =IUCOMP(4HCAVE,IB(JVOLUM+1),NVOLUM)
C      JVO =IB(JVOLUM-IVO)
C      JINL=IB(JVO-6)
C      JINR=IB(JVO-7)
C      IROTL=B(JINL+4)
C      DO 3 I=1,3
C      XL(I)  =B(JINL+I+4)
C      PARL(I)=B(JINL+I+9)
C    3 CONTINUE
C      IROTR=B(JINR+4)
C      DO 4 I=1,3
C      XR(I)  =B(JINR+I+4)
C      PARR(I)=B(JINR+I+9)
C    4 CONTINUE
C      IF(IDEBUG.NE.0)THEN
C        IF(ISWIT(4).NE.0)THEN
C           WRITE(6,*)IROTL,XL(3),PARL(3)
C           WRITE(6,*)IROTR,XR(3),PARR(3)
C        ENDIF
C      ENDIF
C
C         TRANSFORMATION TO LOCAL REFERENCE FRAME
C
C      CALL GFHITS('DRFT','FDOU',1,7,30,0,NUMVS,ITRA,NUMBV,HITS,NHITS)
C
C      IF(NHITS.LE.0)GO TO 99
C      DO 10 IHIT=1,NHITS
C      DO 11 I=1,3
C      X(I) =HITS(I,IHIT)
C      XP(I)=HITS(I+3,IHIT)
C   11 CONTINUE
C      IF(NUMBV(1,IHIT).EQ.1)THEN
C            CALL GITRAN(X,XL,IROTL,XNEW)
C            CALL GITRAN(XP,XD,IROTL,XPNEW)
C            END IF
C      IF(NUMBV(1,IHIT).EQ.2)THEN
C            CALL GITRAN(X,XR,IROTR,XNEW)
C            CALL GITRAN(XP,XD,IROTR,XPNEW)
C            END IF
C      DO 12 I=1,3
C      HITS(I,IHIT)  =XNEW(I)
C      HITS(I+3,IHIT)=XPNEW(I)
C   12 CONTINUE
C   10 CONTINUE
C
C         INTERSECTION WITH SENSE PLANES
C
C      DZPL=2.*PARL(3)/FLOAT(NPL+1)
C      ZPLMIN=-PARL(3)+DZPL
C      DO 20 IHIT=1,NHITS
C      IF(ITRA(IHIT).NE.ITRA(IHIT+1))GO TO 20
C ********  BELOW CHANGE OF COORDINATES TO USE GIPLAN
C      DO 21 I=1,6,3
C      X1(I)  =-HITS(I,IHIT)
C      X1(I+1)= HITS(I+2,IHIT)
C      X1(I+2)= HITS(I+1,IHIT)
C      X2(I)  =-HITS(I,IHIT+1)
C      X2(I+1)= HITS(I+2,IHIT+1)
C      X2(I+2)= HITS(I+1,IHIT+1)
C   21 CONTINUE
C ********
C      ZPL=ZPLMIN
C      DO 22 IPL=1,NPL
C      ZMIN=AMIN1(X1(2),X2(2))
C      ZMAX=AMAX1(X1(2),X2(2))
C      IF((ZPL.LT.ZMIN).OR.(ZPL.GT.ZMAX))GO TO 22
C      CALL GIPLAN(ZPL,X1,X2,S1,S2,1,XINT,SINT,PZINT,IFLAG)
C ********  BACK TO LOCAL REFERENCE FRAME
C ********  AND SAVE INTERACTION COORDINATES
C      XSAV(NUMBV(1,IHIT),IPL,IHIT)=-XINT(1)
C      YSAV(NUMBV(1,IHIT),IPL,IHIT)= XINT(3)
C ********
C      ZPL=ZPL+DZPL
C   22 CONTINUE
C   20 CONTINUE
C
C           DIGITISATION
C
C
   99 RETURN
      END


      SUBROUTINE DLSOUT
C
#include "gcflag.inc"
C
      DIMENSION NUMVS(2),ITRA(50)
      DIMENSION NUMV1(2,50),NUMV2(2,50),NUMV3(2,50),NUMV4(2,50)
      DIMENSION HITS1(9,50),HITS2(9,50),HITS3(9,50),HITS4(9,50)
      DIMENSION VERT(3),PVERT(4)
      SAVE NUMVS
      DATA NUMVS/0,0/
C
      IF(ISWIT(10).LT.3)GO TO 99
      LUN=ISWIT(10)

      CALL GFHITS('DRFT','FSP ',2,9,50,0,NUMVS,ITRA,NUMV1,HITS1,NHITS1)
      CALL GFHITS('DRFT','RSP1',2,9,50,0,NUMVS,ITRA,NUMV2,HITS2,NHITS2)
      CALL GFHITS('DRFT','RSP2',2,9,50,0,NUMVS,ITRA,NUMV3,HITS3,NHITS3)
      CALL GFHITS('DRFT','RSP3',2,9,50,0,NUMVS,ITRA,NUMV4,HITS4,NHITS4)
C
      CALL GFKINE(1,VERT,PVERT,IP,NVERT,UBUF,NUBUF)
C
      WRITE(LUN)NHITS1,NHITS2,NHITS3,NHITS4,IP,VERT,PVERT
      IF(NHITS1.GT.0)THEN
        WRITE(LUN)(NUMV1(2,I),(HITS1(J,I),J=1,9),I=1,NHITS1)
      ENDIF
      IF(NHITS2.GT.0)THEN
        WRITE(LUN)(NUMV2(2,I),(HITS2(J,I),J=1,9),I=1,NHITS2)
      ENDIF
      IF(NHITS3.GT.0)THEN
        WRITE(LUN)(NUMV3(2,I),(HITS3(J,I),J=1,9),I=1,NHITS3)
      ENDIF
      IF(NHITS4.GT.0)THEN
        WRITE(LUN)(NUMV4(2,I),(HITS4(J,I),J=1,9),I=1,NHITS4)
      ENDIF
C
  99  RETURN
      END
