      SUBROUTINE GUKINE
C.
C.    *
C.    *       Generates Kinematics for primary track
C.    *
C.
#include "gcbank.inc"
#include "gcflag.inc"
#include "gckine.inc"
#include "gcunit.inc"
#include "gcscan.inc"
#include "gconsp.inc"
*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEEP,CELOSS.
      COMMON/CELOSS/SEL1(40),SEL1C(40),SER1(40),SER1C(40),SNPAT1(40,4),
     *              SEL2(40),SEL2C(40),SER2(40),SER2C(40),SNPAT2(40,4),
     *              EINTOT,DEDL(40),DEDR(40),FNPAT(40,4)
*KEND.
C

      DIMENSION VERTEX(3),PLAB(3),RNDM(2)
      SAVE VERTEX,PLAB
      DATA VERTEX/3*0./
      DATA PLAB  /3*0./
C.
C.
      IF(SCANFL) THEN
         CALL GSCANK
      ELSE
         CALL VZERO(DEDL,240)
         VERTEX(3) = - 10.*X0 + 0.01
         IF(IKINE.GT.100)THEN
            IK=IKINE-100
            THETA=PKINE(2)*DEGRAD
            PHI=PKINE(3)*DEGRAD
         ELSE
            IK=IKINE
            CALL GRNDM(RNDM,2)
            THETA=PI*RNDM(1)
            PHI=TWOPI*RNDM(2)
         ENDIF
C
         PLAB(1) = PKINE(1)*SIN(THETA)*COS(PHI)
         PLAB(2) = PKINE(1)*SIN(THETA)*SIN(PHI)
         PLAB(3) = PKINE(1)*COS(THETA)
C
         CALL GSVERT(VERTEX,0,0,0,0,NVERT)
         CALL GSKINE(PLAB,IK,NVERT,0,0,NT)

         JK = LQ(JKINE-NT)
         EINTOT = EINTOT + Q(JK+4)


C
C ***          Kinematics debug (controlled by ISWIT(1) )
C
         IF(IDEBUG.EQ.1) THEN
            IF(ISWIT(1).EQ.1) THEN
               CALL GPRINT('VERT',0)
               CALL GPRINT('KINE',0)
            ENDIF
         ENDIF
      ENDIF

C
      END


      SUBROUTINE GUOUT
C.
C.    *
C.    *       User routine called at the end of each event
C.    *
C.
C

#include "gcflag.inc"  
#include "gcscan.inc"  
*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEEP,CELOSS.
      COMMON/CELOSS/SEL1(40),SEL1C(40),SER1(40),SER1C(40),SNPAT1(40,4),
     *              SEL2(40),SEL2C(40),SER2(40),SER2C(40),SNPAT2(40,4),
     *              EINTOT,DEDL(40),DEDR(40),FNPAT(40,4)
*KEND.
      SAVE NID1
      DATA NID1/0/

      IF(ISWIT(9).NE.0) RETURN
      IF(SCANFL) THEN
         CALL GSCANO
      ELSE
         DLC = 0.
         DRC = 0.

         DO 10 I = 1,NL
            SEL1 (I) = SEL1 (I) + DEDL(I)
            SEL2 (I) = SEL2 (I) + DEDL(I)**2

            DLC = DLC + DEDL(I)
            SEL1C(I) = SEL1C(I) + DLC
            SEL2C(I) = SEL2C(I) + DLC**2
   10    CONTINUE

         DO 20 I = 1,NR
            SER1 (I) = SER1 (I) + DEDR(I)
            SER2 (I) = SER2 (I) + DEDR(I)**2

            DRC = DRC + DEDR(I)
            SER1C(I) = SER1C(I) + DRC
            SER2C(I) = SER2C(I) + DRC**2
   20    CONTINUE

         DO 30 IPAT = 1,3
            DO 30 NPL = 1,NL
               SNPAT1(NPL,IPAT) = SNPAT1(NPL,IPAT) + FNPAT(NPL,IPAT)
               SNPAT2(NPL,IPAT) = SNPAT2(NPL,IPAT) + FNPAT(NPL,IPAT)**
     +         2
   30    CONTINUE

         ETOT = 100.*IEVENT*DLC/EINTOT
         CALL HFF1( 1,NID1, ETOT,1.)

      END IF

      END


      SUBROUTINE GUSTEP
C.
C.    *
C.    *       User routine called at the end of each tracking step
C.    *       INWVOL is different from 0 when the track has reached
C.    *              a volume boundary
C.    *       ISTOP is different from 0 if the track has stopped
C.    *
C.
#include "gcbank.inc"  
#include "gctmed.inc"  
#include "gckine.inc"  
#include "gcking.inc"  
#include "gcflag.inc"  
#include "gctrak.inc"  
#include "gcvolu.inc"  
#include "gcscan.inc"  
#include "gconsp.inc"  
#include "gccuts.inc"  
*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEEP,CELOSS.
      COMMON/CELOSS/SEL1(40),SEL1C(40),SER1(40),SER1C(40),SNPAT1(40,4),
     *              SEL2(40),SEL2C(40),SER2(40),SER2C(40),SNPAT2(40,4),
     *              EINTOT,DEDL(40),DEDR(40),FNPAT(40,4)
*KEND.
      DIMENSION NID(100)
      SAVE NID
      DATA NID/100*0/
C
      DATA NROLD / 1 /

c      if (iswit(6).ne.0) call grecord

      IF(SCANFL) THEN
         CALL GSCANU
      ELSE
C
C.
C              Something generated ?

         IF(NGKINE.GT.0) THEN
            DO 10 I=1,NGKINE
               ITYPA = GKIN(5,I)
               CALL GSKING(I)
               IF(ISWIT(9).NE.0)GO TO 10
               ID1=50+ITYPA
               CALL HFF1(ID1,NID(ID1),GKIN(4,I),1.)
               ID2=60+ITYPA
               CALL HFF1(ID2,NID(ID2),GKIN(4,I),1.)
   10       CONTINUE
         ENDIF
         IF(ISWIT(9).NE.0)GO TO 50
         IF(NUMED.EQ.2) THEN
            NRIN = NUMBER(NLEVEL)

C ***          Energy deposited
C
            IF(DESTEP.NE.0.)THEN
               NTUB = NUMBER(NLEVEL-1)
               DEDL(NRIN) = DEDL(NRIN) + DESTEP
               DEDR(NTUB) = DEDR(NTUB) + DESTEP
            ENDIF
C ***          Particle's flux
            IF(SLENG.LE.0.) NROLD = NRIN
            IF(NRIN.NE.NROLD) THEN
               NPL = (NRIN + NROLD)/2 + 1
               FNPAT(NPL,IPART) = FNPAT(NPL,IPART) + 1.
               NROLD = NRIN
            ENDIF


C ***          Process calls frequency
C
            IF(NMEC.LE.0) GO TO 30
            DO 20 IM = 1,NMEC
               IF(LMEC(IM).GT.12) GO TO 20
               IDM = 30 + LMEC(IM)
               CALL HFF1(IDM,NID(IDM),GEKIN,1.)
   20       CONTINUE
   30       CONTINUE

C
C            If track is a photon extrapolate the track
C            directly to the next photon interaction
C            but still computing flux correctly
C
         ENDIF


C             Plot total track length

         IF(ISTOP.NE.0) THEN
            SLENRL = SLENG/X0
            ID3=20+IPART
            CALL HFF1(ID3,NID(ID3),SLENRL,1.)
         ENDIF



C             Debug event

  50     CALL GDEBUG
      END IF

      END


      SUBROUTINE GUTREV
C.
C.    *
C.    *       User routine to control tracking of one event
C.    *
C.    *       Called by GRUN
C.    *
C.
      CALL GTREVE
C
      END


      SUBROUTINE UFILES
*
*            To open FFREAD and HBOOK files
*
      CHARACTER*(*) FILNAM, FSTAT
      PARAMETER (FILNAM='gexam1.dat')
*
      PARAMETER (FSTAT='OLD')
*
      OPEN(UNIT=4,FILE=FILNAM,STATUS=FSTAT,
     +     FORM='FORMATTED')
      END


      SUBROUTINE UGEOM
C
C ***          Define user geometry set up
C

#include "gcbank.inc"
#include "gckine.inc"
*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEND.

      DIMENSION PAR( 8)
      DIMENSION ZLG(6),ALG(6),WLG(6)
      DIMENSION A(3),Z(3),WMAT(3)
      DIMENSION AF(3),ZF(3),WMATF(3)
C
C             Lead glass mixture parameters
C
      DATA ZLG/  82.00,  19.00,  14.00,  11.00,  8.00,  33.00/
      DATA ALG/ 207.19,  39.102,  28.088,  22.99, 15.999,  74.922/
      DATA WLG/ .65994, .00799, .126676, .0040073,.199281, .00200485/
C
C             BGO compound parameters
C
      DATA A/208.98,72.59,15.999/
      DATA Z/83.,32.,8./
      DATA WMAT/4.,3.,12./
C
C             Iron+Nickel+Crome  compound parameters
C
      DATA AF/55.847,58.71,51.998/
      DATA ZF/26.,28.,24./
      DATA WMATF/0.703964,0.099,0.197/
C
C ***          Defines USER particular materials
C
      CALL GSMATE( 9,'ALUMINIUM$', 26.98,13.,2.7   , 8.9,37.2,0,0)
      CALL GSMATE(11,'COPPER$   ', 63.54,29.,8.96  ,1.43,14.8,0,0)
      CALL GSMATE(13,'LEAD$     ',207.19,82.,11.35 ,0.56,18.5,0,0)
      CALL GSMATE(14,'URANIUM$  ',238.03,92.,18.95 ,0.32,12. ,0,0)
      CALL GSMATE(15,'AIR$      ', 14.61,7.3,0.001205,30423.,6750.,0,0)
      CALL GSMATE(16,'VACUUM$ ',1.E-16,1.E-16,1.E-16,1.E+16,1.E+16,0,0)
      CALL GSMIXT(10,'IRON(COMPOUND)$',AF,ZF,7.8,3,WMATF)
      CALL GSMIXT(21,'BGO(COMPOUND)$',A,Z,7.1,-3,WMAT)
      CALL GSMIXT(22,'LEAD GLASS$',ALG,ZLG,5.2,6,WLG)
C
C ***          Defines USER tracking media parameters
C
C

      FIELDM =  0.
      IFIELD =  0
      TMAXFD =  10.
      DMAXMS =  1.
      DEEMAX =  0.05
      EPSIL  =  0.001
      STMIN  =  0.001

      CALL GSTMED( 1,'DEFAULT MEDIUM AIR$'    , 15 , 0 , IFIELD,
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )


      CALL GSTMED( 2,'ABSORBER$'              ,IMAT, 0 , IFIELD,
     *                FIELDM,TMAXFD,DMAXMS,DEEMAX, EPSIL, STMIN, 0 , 0 )


C
C ***          Defines USER'S VOLUMES
C
      NMED1  = 1
      NMED2  = 2
      JMA = LQ(JMATE-IMAT)
      X0 = Q(JMA + 9)
      XR     =  X0/4.

      R1     =  20.   * XR
      R2     =  21.   * XR

      Z1     =  10.    * X0
      Z2     =  11.    * X0


      PAR(1) = 0.
      PAR(2) = R2
      PAR(3) = Z2
      CALL GSVOLU( 'ECAL' , 'TUBE' ,NMED1, PAR , 3 , IVOL )

      PAR(1) = 0.
      PAR(2) = R2
      PAR(3) = 0.5 * X0
      CALL GSVOLU( 'LEAK' , 'TUBE' ,NMED1, PAR , 3 , IVOL )

      PAR(1) = R1
      PAR(2) = R2
      PAR(3) = Z1
      CALL GSVOLU( 'LATR' , 'TUBE' ,NMED1, PAR , 3 , IVOL )

      PAR(1) = 0.
      PAR(2) = R1
      PAR(3) = Z1
      CALL GSVOLU( 'BLOC' , 'TUBE' ,NMED2, PAR , 3 , IVOL )

C
C ***          Position volumes within ECAL
C

      ZC = 0.5 * (Z1 + Z2)
      CALL GSPOS( 'LEAK' ,1, 'ECAL' , 0.   , 0.   , -ZC  , 0,'ONLY')
      CALL GSPOS( 'LEAK' ,2, 'ECAL' , 0.   , 0.   , +ZC  , 0,'ONLY')
      CALL GSPOS( 'LATR' ,1, 'ECAL' , 0.   , 0.   , 0.   , 0,'ONLY')
      CALL GSPOS( 'BLOC' ,1, 'ECAL' , 0.   , 0.   , 0.   , 0,'ONLY')

      CALL GSDVN( 'RTUB' , 'BLOC' ,   NR , 1)
      CALL GSDVN( 'RING' , 'RTUB' ,   NL , 3)


C ***          Close geometry banks. (obligatory system routine)

      CALL GGCLOS

      END


      SUBROUTINE UGINIT
C.
C.    *
C.    *
C.    *        To initialise GEANT/USER  program and read data cards
C.    *
C.
#include "gcunit.inc"
#include "gckine.inc"
#include "gclist.inc"
*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEEP,CELOSS.
      COMMON/CELOSS/SEL1(40),SEL1C(40),SER1(40),SER1C(40),SNPAT1(40,4),
     *              SEL2(40),SEL2C(40),SER2(40),SER2C(40),SNPAT2(40,4),
     *              EINTOT,DEDL(40),DEDR(40),FNPAT(40,4)
*KEND.

C.
*
*             Open user files
*
      CALL UFILES
C.
C             Initialise GEANT
C
      CALL GINIT
      IMAT=10
      PKINE(1)=10.
      PKINE(2)=1.
      EINTOT=0.
      CALL VZERO(SEL1,640)
      NL=20
      NR=20
      CALL FFKEY('BINS',NL,2,'INTEGER')
      CALL FFKEY('MATE',IMAT,1,'INTEGER')

      write(lout,1000)
 1000 format(/,' ========> Reading ffread data cards : type <======='
     +,/,'read 4'
     +,/,'your own data cards if any'
     +,/,'stop',/,'      Now waiting for input',/)

      CALL GFFGO
      CALL GZINIT
      CALL GPART
C
C             Prints version number
C
      WRITE(LOUT,10000)
C

C              Geometry and materials description

      CALL UGEOM
C
      CALL GLOOK('MATE',LPRIN,NPRIN,IM)
      CALL GLOOK('TMED',LPRIN,NPRIN,IT)
      CALL GLOOK('VOLU',LPRIN,NPRIN,IV)
      IF(IM.NE.0)CALL GPRINT('MATE',0)
      IF(IT.NE.0)CALL GPRINT('TMED',0)
      IF(IV.NE.0)CALL GPRINT('VOLU',0)

C              Energy loss and cross-sections initialisations

      CALL GPHYSI

C             Define user histograms

      CALL UHINIT

10000 FORMAT(/,'  GEXAM1 VERSION 1.00 ',/)

      END


      SUBROUTINE UGLAST
C.    *
C.    *
C.    *      Termination routine to print histograms and statistics
C.    *
C.    *

#include "gcflag.inc"  
*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEEP,CELOSS.
      COMMON/CELOSS/SEL1(40),SEL1C(40),SER1(40),SER1C(40),SNPAT1(40,4),
     *              SEL2(40),SEL2C(40),SER2(40),SER2C(40),SNPAT2(40,4),
     *              EINTOT,DEDL(40),DEDR(40),FNPAT(40,4)
*KEND.
      COMMON/SCLAST/XSEL1(40),XSEL2(40),XSEL1C(40),XSEL2C(40)
     +             ,XSER1(40),XSER2(40),XSER1C(40),XSER2C(40)

      DIMENSION PAT1(50),PAT2(50)
C
      IF(ISWIT(10).EQ.0)CALL GLAST

C ***          Normalize and print energy distribution


      if (eintot.eq.0) return;
      CALL VZERO(XSEL1,320)
      CALL VZERO(PAT1,50)
      CALL VZERO(PAT2,50)
      CNORM  = 100./EINTOT
      XEVENT=IEVENT
C
      DO 10 I = 1,NL
      XSEL1 (I) = CNORM * SEL1 (I)
      XSEL2 (I) = CNORM*SQRT(ABS(XEVENT*SEL2 (I) - SEL1 (I)**2))

      XSEL1C(I) = CNORM * SEL1C(I)
      XSEL2C(I) = CNORM*SQRT(ABS(XEVENT*SEL2C(I) - SEL1C(I)**2))
   10 CONTINUE

      CALL HPAK (2,XSEL1 )
      CALL HPAKE(2,XSEL2 )

      CALL HPAK (4,XSEL1C)
      CALL HPAKE(4,XSEL2C)

      DO 20 I = 1,NR
      XSER1 (I) = CNORM * SER1 (I)
      XSER2 (I) = CNORM*SQRT(ABS(XEVENT*SER2 (I) - SER1 (I)**2))

      XSER1C(I) = CNORM * SER1C(I)
      XSER2C(I) = CNORM*SQRT(ABS(XEVENT*SER2C(I) - SER1C(I)**2))
   20 CONTINUE

      CALL HPAK (3,XSER1 )
      CALL HPAKE(3,XSER2 )

      CALL HPAK (5,XSER1C)
      CALL HPAKE(5,XSER2C)

      DO 40 IP = 1,3
      DO 30 I  = 1,NL
      PAT1(I) = SNPAT1(I,IP) / XEVENT
      PAT2(I) = SQRT(SNPAT2(I,IP)/XEVENT - PAT1(I)**2)
   30 CONTINUE
      CALL HPAK (10+IP,PAT1)
      CALL HPAKE(10+IP,PAT2)
   40 CONTINUE

      PRINT 10000
      PRINT 10200,( XSEL2 (I),I=1,NL)
      PRINT 10300,( XSEL2C(I),I=1,NL)
      PRINT 10100
      PRINT 10200,( XSER2 (I),I=1,NR)
      PRINT 10300,( XSER2C(I),I=1,NR)
10000 FORMAT(///,40X,'LONGITUDINAL PROFIL',/)
10100 FORMAT(///,40X,'   RADIAL PROFIL   ',/)
10200 FORMAT(//,30X,'ERROR ON PROFIL     VALUES',/
     +         ,(10X,10F10.4))
10300 FORMAT(//,30X,'ERROR ON CUMULATIVE VALUES',/
     +         ,(10X,10F10.4))

C
C             Save histograms
C
      CALL HRPUT(0,'gexam1.hist',' ')
C
      IF(ISWIT(10).EQ.0)THEN
         CALL HIDOPT(0,'BLAC')
         CALL HISTDO
      ENDIF
      END


      SUBROUTINE UHINIT
C
C     *       To book the user's histograms
C     *

*KEEP,PVOLUM.
      COMMON/PVOLUM/ NL,NR,IMAT,X0
*KEND.

C
C ***            Histograms for shower development
C
      NBINZ=NL+1
      NBINR=NR+1
      ZMAX=NBINZ+1
      RMAX=NBINR+1
C
      CALL HBOOK1(1,'TOTAL ENERGY DEPOSITION$'
     *,100, 81.,101., 0.0)
      CALL HBOOK1(2,'LONGIT ENERGY DEPOSITION $'
     *, NBINZ, 1.,ZMAX,  0.0)
      CALL HBOOK1(4,'CUMUL LONGIT ENERGY DEP. $'
     *, NBINZ, 1.,ZMAX,  0.0)
      CALL HBOOK1(3,'RADIAL ENERGY DEPOSITION $'
     *, NBINR, 1.,RMAX,  0.0)
      CALL HBOOK1(5,'CUMUL RADIAL ENERGY DEP. $'
     *, NBINR, 1.,RMAX,  0.0)
      CALL HBOOK1(11,'NB OF GAMMA PER PLANE$'
     *, NBINZ, 1.,ZMAX,  0.0)
      CALL HBOOK1(12,'NB OF E +   PER PLANE$'
     *, NBINZ, 1.,ZMAX,  0.0)
      CALL HBOOK1(13,'NB OF E -   PER PLANE$'
     *, NBINZ, 1.,ZMAX,  0.0)

      CALL HBIGBI(0,4)


C
C ***          Histograms for detailed studies
C
      CALL HBOOK1(21,'TOTAL GAMMA LENGHT IN RL$'
     *,100, 0.,  10. , 0.)
      CALL HBOOK1(22,'TOTAL POSIT LENGHT IN RL$'
     *,100, 0., 5.   , 0.)
      CALL HBOOK1(23,'TOTAL ELECT LENGHT IN RL$'
     *,100, 0., 5.   , 0.)

      CALL HBOOK1(31,'NUMBER OF NEXT CALLS$'
     *,100, 0., 0.1, 0.)
      CALL HBOOK1(32,'NUMBER OF MULS CALLS$'
     *,100, 0., 0.01, 0.)
      CALL HBOOK1(33,'NUMBER OF LOSS CALLS$'
     *,100, 0., 0.01, 0.)
      CALL HBOOK1(34,'NUMBER OF FIEL CALLS$'
     *,100, 0., 0.1, 0.)
      CALL HBOOK1(35,'NUMBER OF DCAY CALLS$'
     *,100, 0., 0.1, 0.)
      CALL HBOOK1(36,'NUMBER OF PAIR CALLS$'
     *,100, 0., 0.1, 0.)
      CALL HBOOK1(37,'NUMBER OF COMP CALLS$'
     *,100, 0., 0.01, 0.)
      CALL HBOOK1(38,'NUMBER OF PHOT CALLS$'
     *,100, 0., 0.01, 0.)
      CALL HBOOK1(39,'NUMBER OF BREM CALLS$'
     *,100, 0., 0.1, 0.)
      CALL HBOOK1(40,'NUMBER OF DRAY CALLS$'
     *,100, 0., 0.1, 0.)
      CALL HBOOK1(41,'NUMBER OF ANNI CALLS$'
     *,100, 0., 0.01, 0.)

      CALL HBOOK1(51,'ENERGY DISTR OF GAMMAS$'
     *,100,0.,0.1,0.)
      CALL HBOOK1(52,'ENERGY DISTR OF POSITRONS$'
     *,100,0.,0.1,0.)
      CALL HBOOK1(53,'ENERGY DISTR OF ELECTRONS$'
     *,100,0.,0.1,0.)
      CALL HBOOK1(61,'ENERGY DISTR OF GAMMAS$'
     *,100,0.,0.01,0.)
      CALL HBOOK1(62,'ENERGY DISTR OF POSITRONS$'
     *,100,0.,0.01,0.)
      CALL HBOOK1(63,'ENERGY DISTR OF ELECTRONS$'
     *,100,0.,0.01,0.)



      END
