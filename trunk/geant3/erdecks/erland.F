*
* $Id: erland.F 229 2008-11-17 18:13:31Z brun $
*
* $Log: erland.F,v $
* Revision 1.2  2007/03/31 13:47:09  brun
*  From Andrea Fontana and Alberto Rotondi
* Added new calculation for sigma of straggling in energy loss to include
* in Geane the Urban/Landau approximation, as explained in the Geant
* manual and related papers.
* The model parametrization can be controlled with a new user variable
* (GCALPHA) in GCCUTS common block: 1 is for old gaussian model valid for
* dense materials, other values (see the enclosed report) are for gaseous
* materials.
*
* Revision 1.1.1.1  2002/07/24 15:56:26  rdm
* initial import into CVS
*
* Revision 1.1.1.1  2002/06/16 15:18:35  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/03/06 15:37:34  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERLAND (STEP,ZMAT,AMAT,RHO,P,E,XMASS,DEDX,DEDX2)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates energy straggling using Gaussian theory in a step  *
C.    *                                                                *
C.    *  Input  STEP   =  current step-length (cm)                     *
C.    *         ZMAT   =  Z of the material                            *
C.    *         AMAT   =  A of the material                            *
C.    *         RHO    =  density of the material                      *
C.    *         P      =  momentum of the particle                     *
C.    *         E      =  energy   of the particle                     *
C.    *         XMASS  =  mass     of the particle                     *
C.    *         DEDX   =  mean energy loss in STEP                     *
C.    *                                                                *
C.    *  Output DEDX2  =  mean square of the straggling in G           *
C.    *                                                                *
C.    *    ==>Called by : ERTRCH                                       *
C.    *       Author    E.Nagy  *********                              *
C.    *       Original routine : GLANDO                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconst.inc"
c included gcmate COMMON to get radiation length
#include "geant321/gcmate.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcmore.inc"
*
      REAL*4 CMED, DE, ION, RU, FU1, FU2, EU1, EU2, EU2F, DEDX
      REAL*4 EUMED,SIG1, SIG2, SIG3, SIG4,CN1, CN2, CN3, XUMASS 
      REAL*4 EUMED2, E99, A99 
      REAL*4 SALPHA, RLAMAX, RLAMED
      REAL*4 LX,EMB,S2B,SB,E1B,E2B,DEDXB
*
      IF (STEP.LT.1.E-7.OR.P.LT.1.E-10) THEN
         DEDX2=0.
         RETURN
      ENDIF
*
*     Calculate xi factor (KeV).
      BETA   = P/E
      GAMMA  = E/XMASS
      XI     = (153.5*ZMAT*STEP*RHO)/(AMAT*BETA*BETA)
*
*     Maximum energy transfer to atomic electron (KeV).
      ETA    = BETA*GAMMA
      ETASQ  = ETA*ETA
      RATIO  = EMASS/XMASS
      F1     = 2.*EMASS*ETASQ
      F2     = 1.+2.*RATIO*GAMMA+RATIO*RATIO
      EMAX   = 1.E+6*F1/F2
       
*     old GEANE standard  straggling sigma in Gev**2
      A99 = GCALPHA    
      IF (A99.EQ.1.) THEN  
          DEDX2 = XI*EMAX*(1.-(BETA*BETA/2.))*1.E-12
          RETURN           	   		   
      ENDIF
      
C      new calculation for sigma of straggling
C      by A. Fontana and A. Rotondi  March 2007
C      revised in september 2008
C
      SALPHA=0.
      IF(XI/EMAX .GT. 0.01)  THEN
C           Vavilov-Gaussian regime:  sigma^2(E)  in GeV^2
         DEDX2 = XI*EMAX*(1.-(BETA*BETA/2.))*1.E-12
      ELSE 
C
C          Urban/Landau approximation
C
C
C          energies in eV
C
         ION=16.*ZMAT**0.9
         RU=0.6
         FU2= 0
         IF(ZMAT .GT. 2) FU2=2./ZMAT
         FU1= 1.-FU2
         EU2 =10.*ZMAT*ZMAT
         EU2F=EU2**FU2
         EU1= (ION/EU2F)**(1./FU1)  

C          energies in keV
         ION=ION/1000.
         EU1=EU1/1000.
         EU2=EU2/1000.
         XUMASS = XMASS*1.E+06
         CMED=DEDX*1.E+06/STEP
	 E99 = ION/(1.- A99*EMAX/(EMAX+ION))
         
         SIG1=(CMED*FU1/EU1)*(ALOG(2.*XUMASS*((BETA*GAMMA)**2)/EU1)
     +        - BETA*BETA)
         SIG1= (1-RU)*SIG1/(ALOG(2.*XUMASS*((BETA*GAMMA)**2)/ION)
     +        -BETA*BETA)
         SIG2=(CMED*FU2/EU2)*(ALOG(2.*XUMASS*((BETA*GAMMA)**2)/EU2) 
     +        - BETA*BETA)
         SIG2= (1-RU)*SIG2/(ALOG(2.*XUMASS*((BETA*GAMMA)**2)/ION)
     +        -BETA*BETA)
         SIG3= CMED*RU*EMAX/(ION*(EMAX+ION)*ALOG((EMAX+ION)/ION))
C         number of collisions
         CN1=SIG1*STEP
         CN2=SIG2*STEP
         CN3=SIG3*STEP
         IF(CN1+CN2+CN3 .LT. 50.) THEN
C            Urban model
            EUMED= (ION*(EMAX+ION)/EMAX)*ALOG(E99/ION)
	    EUMED2 = (ION*(EMAX+ION)/EMAX)*(E99-ION)
	    SIG4= EUMED2 - EUMED*EUMED
            DEDX2=CN1*EU1*EU1 + CN2*EU2*EU2 + CN3*EUMED*EUMED 
     +             + CN3*SIG4	     
         ELSE
C
C            truncated Landau distribution 
C              see GEANT3 manual W5013
            RLAMED = -0.422784 -BETA*BETA - ALOG(XI/EMAX)
            RLAMAX =  0.60715 + 1.1934*RLAMED 
     +     +(0.67794 + 0.052382*RLAMED)*EXP(0.94753+0.74442*RLAMED)

C            from lambda max to SALPHA=sigma (empirical polynomial)
            IF(RLAMAX .LE. 1010) THEN
               SALPHA =    1.975560  
     +              +9.898841e-02 *RLAMAX  
     +              -2.828670e-04 *RLAMAX**2  
     +              +5.345406e-07 *RLAMAX**3   
     +              -4.942035e-10 *RLAMAX**4  
     +              +1.729807e-13 *RLAMAX**5 

            ELSE 
               SALPHA =1.871887E+01 + 1.296254E-02 *RLAMAX  
            END IF
            
C           alpha=54.6  corresponds to a 0.9996 maximum cut
            IF(SALPHA > 54.6) SALPHA=54.6  
            DEDX2 = (SALPHA*XI)**2         
         END IF
C         sigma^2(E) in  GeV^2
         DEDX2 = DEDX2*1.E-12
      END IF
C      
C  inclusion of bremsstrahlung fluctuations DEDXB for electrons
C  22.09.2008 by Alberto Rotondi
C      
      IF(XMASS.LT.0.6E-3) THEN
         LX = 1.442695*STEP/RADL
         EMB=E/(2**LX)
	 S2B=E*E*(1/3**LX - 1/4**LX)
	 SB=SQRT(ABS(S2B))
         DEDXB = 1.2*SB
         DEDX2 = DEDX2 + DEDXB**2
      ENDIF
C      PRINT *,XI/EMAX,AMAT,ZMAT,RHO
c      PRINT *, SALPHA, A99, CN1+CN2+CN3, SQRT(DEDX2)
*
      END
