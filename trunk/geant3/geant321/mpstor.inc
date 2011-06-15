*
* $Id: mpstor.inc 220 2007-11-19 16:08:06Z rdm $
*
* $Log: mpstor.inc,v $
* Revision 1.1.1.1  2002/06/16 15:18:38  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:17  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:20:45  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_MPSTOR_INC
#define CERNLIB_GEANT321_MPSTOR_INC
*
*
* mpstor.inc
*
      PARAMETER(IDNEU  = 1)
      PARAMETER(IDHEVY = 2)
      PARAMETER(IDGAMA = 3)
      COMMON/ MPSTOR / EN(MAXPAR),UN(MAXPAR),VN(MAXPAR),WN(MAXPAR),
     +                 AGEN(MAXPAR),MTN(MAXPAR),AMN(MAXPAR),
     +                 ZMN(MAXPAR),IDN(MAXPAR),
     +                 EP,UP,VP,WP,MTP,AGEP,AMP,ZMP,
     +                 NNEU,NHEVY,NGAMA,NPSTOR

#endif
