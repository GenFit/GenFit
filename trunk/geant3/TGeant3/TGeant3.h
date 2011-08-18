#ifndef ROOT_TGeant3
#define ROOT_TGeant3
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: TGeant3.h 247 2009-10-09 13:20:34Z brun $ */

////////////////////////////////////////////////
//  C++ interface to Geant3 basic routines    //
////////////////////////////////////////////////

#define WITHG3
#ifdef WITHROOT
#undef WITHG3
#endif
#ifdef WITHBOTH
#undef WITHG3
#undef WITHROOT
#endif

#include "TVirtualMC.h"
#include "TMCProcess.h"
#include "TMCParticleType.h"
#include "TGeoMCGeometry.h"
#include "TObjArray.h"
#include "TArrayI.h"

class TGeoHMatrix;
class TArrayD;
class TString;

//______________________________________________________________
//
//       Geant3 prototypes for commons
//
//______________________________________________________________
//

//----------GCONST
//     COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
typedef struct {
   Float_t   pi; 
   Float_t   twopi;
   Float_t   piby2;
   Float_t   degrad;
   Float_t   raddeg;
   Float_t   clight;
   Float_t   big;
   Float_t   emass;
}Gconst_t;

//----------GCONSX
//      COMMON/GCONSX/EMMU,PMASS,AVO
typedef struct {
   Float_t   emmu;
   Float_t   pmass;
   Float_t   avo;
} Gconsx_t;

//---------- GCJUMP
//     PARAMETER    (MAXJMP=30)
//     COMMON/GCJUMP/JUDCAY, JUDIGI, JUDTIM, JUFLD , JUHADR, JUIGET,
//     +              JUINME, JUINTI, JUKINE, JUNEAR, JUOUT , JUPHAD,
//     +              JUSKIP, JUSTEP, JUSWIM, JUTRAK, JUTREV, JUVIEW,
//     +              JUPARA
//      DIMENSION     JMPADR(MAXJMP)
//      EQUIVALENCE  (JMPADR(1), JUDCAY)
typedef struct {
   Int_t   judcay; 
   Int_t   judigi; 
   Int_t   judtim; 
   Int_t   jufld ; 
   Int_t   juhadr; 
   Int_t   juiget;
   Int_t   juinme; 
   Int_t   juinti; 
   Int_t   jukine; 
   Int_t   junear; 
   Int_t   juout ; 
   Int_t   juphad;
   Int_t   juskip; 
   Int_t   justep; 
   Int_t   juswim; 
   Int_t   jutrak; 
   Int_t   jutrev; 
   Int_t   juview;
   Int_t   jupara;
   Int_t   jmpadr[30];
} Gcjump_t;

//----------QUEST
//      COMMON/QUEST/IQUEST(100)
typedef struct {
  Int_t   iquest[100];
} Quest_t;

//----------GCBANK
//      COMMON/GCBANK/NZEBRA,GVERSN,ZVERSN,IXSTOR,IXDIV,IXCONS,FENDQ(16)
//     +             ,LMAIN,LR1,WS(KWBANK)
typedef struct {
  Int_t nzebra;
  Float_t gversn;
  Float_t zversn;
  Int_t ixstor;
  Int_t ixdiv;
  Int_t ixcons;
  Float_t fendq[16];
  Int_t lmain;
  Int_t lr1;
} Gcbank_t;

//----------GCLINK
//      COMMON/GCLINK/JDIGI ,JDRAW ,JHEAD ,JHITS ,JKINE ,JMATE ,JPART
//     +      ,JROTM ,JRUNG ,JSET  ,JSTAK ,JGSTAT,JTMED ,JTRACK,JVERTX
//     +      ,JVOLUM,JXYZ  ,JGPAR ,JGPAR2,JSKLT
typedef struct {
  Int_t    jdigi;
  Int_t    jdraw;
  Int_t    jhead;
  Int_t    jhits;
  Int_t    jkine;
  Int_t    jmate;
  Int_t    jpart;
  Int_t    jrotm;
  Int_t    jrung;
  Int_t    jset;
  Int_t    jstak;
  Int_t    jgstat;
  Int_t    jtmed;
  Int_t    jtrack;
  Int_t    jvertx;
  Int_t    jvolum;
  Int_t    jxyz;
  Int_t    jgpar;
  Int_t    jgpar2;
  Int_t    jsklt;
} Gclink_t;


//----------GCFLAG
//      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
//     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
typedef struct {
  Int_t    idebug;
  Int_t    idemin;
  Int_t    idemax;
  Int_t    itest;
  Int_t    idrun;
  Int_t    idevt;
  Int_t    ieorun;
  Int_t    ieotri;
  Int_t    ievent;
  Int_t    iswit[10];
  Int_t    ifinit[20];
  Int_t    nevent;
  Int_t    nrndm[2];
} Gcflag_t;

//----------GCKINE
//      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
//     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
typedef struct {
  Int_t    ikine;
  Float_t  pkine[10];
  Int_t    itra;
  Int_t    istak;
  Int_t    ivert;
  Int_t    ipart;
  Int_t    itrtyp;
  Int_t    napart[5];
  Float_t  amass;
  Float_t  charge;
  Float_t  tlife;
  Float_t  vert[3];
  Float_t  pvert[4];
  Int_t    ipaold;
} Gckine_t;

//----------GCKING
//      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN),
//     +                           TOFD(MXGKIN),IFLGK(MXGKIN)
#define MXGKIN 100
typedef struct  {
  Int_t    kcase;
  Int_t    ngkine;
  Float_t  gkin[MXGKIN][5];
  Float_t  tofd[MXGKIN];
  Int_t    iflgk[MXGKIN];
} Gcking_t;

//----------GCKIN2
//      COMMON/GCKIN2/NGPHOT,XPHOT(11,MXPHOT)
#define MXPHOT 800
typedef struct {
  Int_t ngphot;
  Float_t xphot[MXPHOT][11];
} Gckin2_t;

//----------GCKIN3
//      COMMON/GCKIN3/GPOS(3,MXGKIN)
typedef struct {
  Float_t gpos[MXGKIN][3];
} Gckin3_t;

//----------GCMATE
//      COMMON/GCMATE/NMAT,NAMATE(5),A,Z,DENS,RADL,ABSL
typedef struct {
  Int_t    nmat;
  Int_t    namate[5];
  Float_t  a;
  Float_t  z;
  Float_t  dens;
  Float_t  radl;
  Float_t  absl;
} Gcmate_t;

//----------GCTMED
//      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX
//     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD
typedef struct {
  Int_t    numed;
  Int_t    natmed[5];
  Int_t    isvol;
  Int_t    ifield;
  Float_t  fieldm;
  Float_t  tmaxfd;
  Float_t  stemax;
  Float_t  deemax;
  Float_t  epsil;
  Float_t  stmin;
  Float_t  cfield;
  Float_t  prec;
  Int_t    iupd;
  Int_t    istpar;
  Int_t    numold;
} Gctmed_t;

//----------GCTRAK
#define MAXMEC 30
//      PARAMETER (MAXMEC=30)
//      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
//     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
//     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
//     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
//     + ,NLVSAV,ISTORY
typedef struct {
  Float_t  vect[7];
  Float_t  getot;
  Float_t  gekin;
  Float_t  vout[7];
  Int_t    nmec;
  Int_t    lmec[MAXMEC];
  Int_t    namec[MAXMEC];
  Int_t    nstep;
  Int_t    maxnst;
  Float_t  destep;
  Float_t  destel;
  Float_t  safety;
  Float_t  sleng;
  Float_t  step;
  Float_t  snext;
  Float_t  sfield;
  Float_t  tofg;
  Float_t  gekrat;
  Float_t  upwght;
  Int_t    ignext;
  Int_t    inwvol;
  Int_t    istop;
  Int_t    igauto;
  Int_t    iekbin;
  Int_t    ilosl;
  Int_t    imull;
  Int_t    ingoto;
  Int_t    nldown;
  Int_t    nlevin;
  Int_t    nlsav;
  Int_t    istory;
} Gctrak_t;

//----------GCVOLU
//      COMMON/GCVOLU/NLEVEL,NAMES(15),NUMBER(15),
//     +LVOLUM(15),LINDEX(15),INFROM,NLEVMX,NLDEV(15),LINMX(15),
//     +GTRAN(3,15),GRMAT(10,15),GONLY(15),GLX(3)
typedef struct {
  Int_t    nlevel;
  Int_t    names[15];
  Int_t    number[15];
  Int_t    lvolum[15];
  Int_t    lindex[15];
  Int_t    infrom;
  Int_t    nlevmx;
  Int_t    nldev[15];
  Int_t    linmx[15];
  Float_t  gtran[15][3];
  Float_t  grmat[15][10];
  Float_t  gonly[15];
  Float_t  glx[3];
} Gcvolu_t;

//----------GCSETS
//  COMMON/GCSETS/IHSET,IHDET,ISET,IDET,IDTYPE,NVNAME,NUMBV(20)
typedef struct {
  Int_t    ihset;
  Int_t    ihdet;
  Int_t    iset;
  Int_t    idet;
  Int_t    idtype;
  Int_t    nvname;
  Int_t    numbv[20];
} Gcsets_t;

//----------GCNUM
//   COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
//  +            ,NSTMAX,NVERTX,NHEAD,NBIT
typedef struct {
  Int_t    nmate;
  Int_t    nvolum;
  Int_t    nrotm;
  Int_t    ntmed;
  Int_t    ntmult;
  Int_t    ntrack;
  Int_t    npart;
  Int_t    nstmax;
  Int_t    nvertx;
  Int_t    nhead;
  Int_t    nbit;
} Gcnum_t;

//----------GCCUTS
//  COMMON/GCCUTS/CUTGAM,CUTELE,CUTNEU,CUTHAD,CUTMUO,BCUTE,BCUTM
//   +             ,DCUTE ,DCUTM ,PPCUTM,TOFMAX,GCUTS(5)
typedef struct {
  Float_t cutgam;
  Float_t cutele;
  Float_t cutneu;
  Float_t cuthad;
  Float_t cutmuo;
  Float_t bcute;
  Float_t bcutm;
  Float_t dcute;
  Float_t dcutm;
  Float_t ppcutm;
  Float_t tofmax;
  Float_t gcuts[5];
} Gccuts_t;

//----------GCMORE
//      COMMON/GCMORE/GCALPHA,ICLOSE,PFINAL(3),DSTRT,WIRE1(3),WIRE2(3),
//     +              P1(3),P2(3),P3(3),CLENG(3)
typedef struct {
  Float_t  gcalpha;
  Int_t    iclose;
  Float_t  pfinal[3];
  Float_t  dstrt;
  Float_t  wire1[3];
  Float_t  wire2[3];
  Float_t  p1[3];
  Float_t  p2[3];
  Float_t  p3[3];
  Float_t  cleng[3];
} Gcmore_t;

//----------GCMULO
//      COMMON/GCMULO/SINMUL(101),COSMUL(101),SQRMUL(101),OMCMOL,CHCMOL
//     +  ,EKMIN,EKMAX,NEKBIN,NEK1,EKINV,GEKA,GEKB,EKBIN(200),ELOW(200)
typedef struct {
  Float_t sinmul[101];
  Float_t cosmul[101];
  Float_t sqrmul[101];
  Float_t omcmol;
  Float_t chcmol;
  Float_t ekmin;
  Float_t ekmax;
  Int_t   nekbin;
  Int_t   nek1;
  Float_t ekinv;
  Float_t geka;
  Float_t gekb;
  Float_t ekbin[200];
  Float_t elow[200];
} Gcmulo_t;

//----------GCPHYS
//      COMMON/GCPHYS/IPAIR,SPAIR,SLPAIR,ZINTPA,STEPPA
//     +             ,ICOMP,SCOMP,SLCOMP,ZINTCO,STEPCO
//     +             ,IPHOT,SPHOT,SLPHOT,ZINTPH,STEPPH
//     +             ,IPFIS,SPFIS,SLPFIS,ZINTPF,STEPPF
//     +             ,IDRAY,SDRAY,SLDRAY,ZINTDR,STEPDR
//     +             ,IANNI,SANNI,SLANNI,ZINTAN,STEPAN
//     +             ,IBREM,SBREM,SLBREM,ZINTBR,STEPBR
//     +             ,IHADR,SHADR,SLHADR,ZINTHA,STEPHA
//     +             ,IMUNU,SMUNU,SLMUNU,ZINTMU,STEPMU
//     +             ,IDCAY,SDCAY,SLIFE ,SUMLIF,DPHYS1
//     +             ,ILOSS,SLOSS,SOLOSS,STLOSS,DPHYS2
//     +             ,IMULS,SMULS,SOMULS,STMULS,DPHYS3
//     +             ,IRAYL,SRAYL,SLRAYL,ZINTRA,STEPRA
typedef struct {
  Int_t    ipair;
  Float_t  spair;
  Float_t  slpair;
  Float_t  zintpa;
  Float_t  steppa;
  Int_t    icomp;
  Float_t  scomp;
  Float_t  slcomp;
  Float_t  zintco;
  Float_t  stepco;
  Int_t    iphot;
  Float_t  sphot;
  Float_t  slphot;
  Float_t  zintph;
  Float_t  stepph;
  Int_t    ipfis;
  Float_t  spfis;
  Float_t  slpfis;
  Float_t  zintpf;
  Float_t  steppf;
  Int_t    idray;
  Float_t  sdray;
  Float_t  sldray;
  Float_t  zintdr;
  Float_t  stepdr;
  Int_t    ianni;
  Float_t  sanni;
  Float_t  slanni;
  Float_t  zintan;
  Float_t  stepan;
  Int_t    ibrem;
  Float_t  sbrem;
  Float_t  slbrem;
  Float_t  zintbr;
  Float_t  stepbr;
  Int_t    ihadr;
  Float_t  shadr;
  Float_t  slhadr;
  Float_t  zintha;
  Float_t  stepha;
  Int_t    imunu;
  Float_t  smunu;
  Float_t  slmunu;
  Float_t  zintmu;
  Float_t  stepmu;
  Int_t    idcay;
  Float_t  sdcay;
  Float_t  slife;
  Float_t  sumlif;
  Float_t  dphys1;
  Int_t    iloss;
  Float_t  sloss;
  Float_t  soloss;
  Float_t  stloss;
  Float_t  dphys2;
  Int_t    imuls;
  Float_t  smuls;
  Float_t  somuls;
  Float_t  stmuls;
  Float_t  dphys3;
  Int_t    irayl;
  Float_t  srayl;
  Float_t  slrayl;
  Float_t  zintra;
  Float_t  stepra;
} Gcphys_t;

//----------GCPHLT
//      COMMON/GCPHLT/ILABS,SLABS,SLLABS,ZINTLA,STEPLA
//     +             ,ISYNC
//     +             ,ISTRA
typedef struct {
  Int_t ilabs;
  Float_t slabs;
  Float_t sllabs;
  Float_t zintla;
  Float_t stepla;
  Int_t isync;
  Int_t istra;
} Gcphlt_t;

//----------GCOPTI
//      COMMON/GCOPTI/IOPTIM
typedef struct {
  Int_t   ioptim;
} Gcopti_t;

//----------GCTLIT
//      COMMON/GCTLIT/THRIND,PMIN,DP,DNDL,JMIN,ITCKOV,IMCKOV,NPCKOV
typedef struct {
  Float_t   thrind;
  Float_t   pmin;
  Float_t   dp;
  Float_t   dndl;
  Int_t     jmin;
  Int_t     itckov;
  Int_t     imckov;
  Int_t     npckov;
} Gctlit_t;

//----------GCVDMA
//      COMMON/GCVDMA/NVMANY,MANYLE(20),MANYNA(20,15),
//     +MANYNU(20,15),NFMANY,MYCOUN,IMYSE,RAYTRA,VECCOS(3)
typedef struct {
  Int_t     vdma[624];
  Float_t   raytra;
  Float_t   veccos[3];
} Gcvdma_t;

//----------GCTPOL
#define MAXME1 30
//      COMMON/GCTPOL/POLAR(3), NAMEC1(MAXME1)
typedef struct {
  Float_t polar[3];
  Int_t   namec1[MAXME1];
} Gctpol_t;

/************************************************************************
 *                                                                      *
 *      Commons for GEANE                                               *
 *                                                                      *
 ************************************************************************/

//------------ERTRIO
//    INTEGER          MXPRED
//    PARAMETER (MXPRED = 10)
//    DOUBLE PRECISION ERDTRP
//    REAL             ERRIN, ERROUT, ERTRSP, ERXIN, ERXOUT, ERPIN,
//   +                 ERPOUT
//    INTEGER          NEPRED, INLIST, ILPRED, IEPRED
//    COMMON /ERTRIO/  ERDTRP(5,5,MXPRED), ERRIN(15), ERROUT(15,MXPRED),
//   +                 ERTRSP(5,5,MXPRED), ERXIN( 3), ERXOUT( 3,MXPRED),
//   +                 ERPIN(3), ERPOUT(3,MXPRED), NEPRED,INLIST,ILPRED,
//   +                 IEPRED(MXPRED)
//

#define MXPRED 10
typedef struct {
  Double_t erdtrp[MXPRED*5*5];
  Float_t  errin[15];
  Float_t  errout[MXPRED*15];
  Float_t  ertrsp[MXPRED*5*5];
  Float_t  erxin[3];
  Float_t  erxout[MXPRED*3];
  Float_t  erpin[3];
  Float_t  erpout[MXPRED*3];
  Int_t    nepred;
  Int_t    inlist;
  Int_t    ilpred;
  Int_t    iepred;
} Ertrio_t;

typedef struct {
  Int_t iertr;
  Int_t iertr1;
  Int_t iertr2;
} Ertrio1_t;

//-----------EROTPS
//    CHARACTER*8     CHOPTI
//    LOGICAL         LEEXAC, LELENG, LEONLY, LEPLAN, LEPOIN, LEVOLU
//    REAL            ERPLI, ERPLO, ERLENG
//    INTEGER         NAMEER, NUMVER, IOVLER
//    COMMON /EROPTS/ ERPLI(3,2), ERPLO(3,4,MXPRED), ERLENG(MXPRED),
//   +                NAMEER(MXPRED), NUMVER(MXPRED), IOVLER(MXPRED),
//   +                LEEXAC, LELENG, LEONLY, LEPLAN, LEPOIN, LEVOLU
//    COMMON /EROPTC/CHOPTI

typedef struct {
  Float_t   erpli[3*2];
  Float_t   erplo[MXPRED*3*4];
  Float_t   erleng[MXPRED];
  Int_t     nameer[MXPRED];
  Int_t     numver[MXPRED];
  Int_t     iovler[MXPRED];
  Int_t    leexac;
  Int_t    leleng;
  Int_t    leonly;
  Int_t    leplan;
  Int_t    lepoin;
  Int_t    levolu;
} Eropts_t;

typedef struct {
  char chopti[8];
} Eroptc_t;

//-------TRCOM3: A. Panzarasa
// COMMON /TRCOM3/ A(5,5),B(5,5),S(15),TN(3),T(5),
//          COSL,SINL,COSP,SINP,COSLI,NEW
//
typedef struct {
  Double_t a[5][5];
  Double_t b[5][5];
  Double_t s[15];
  Double_t tn[3];
  Double_t t[5];
  Double_t cosl;
  Double_t sinl;
  Double_t cosp;
  Double_t sinp;
  Double_t cosl1;
  Int_t NEW;
} Trcom3_t;


//-------ERWORK
//    DOUBLE PRECISION EI, EF, ASDSC
//    COMMON /ERWORK/ EI(15), EF(15), ASDSC(5,5),
//   +                   XI(3), PPI(3), HI(9),
//   +                   XF(3), PF(3),  HF(9),
//   +                   CHTR, DEDX2, BACKTR, CUTEK, TLGCM2, TLRAD

typedef struct {
  Double_t  ei[15];
  Double_t  ef[15];
  Double_t  asdsc[5*5];
  Float_t   xi[3];
  Float_t   ppi[3];
  Float_t   hi[9];
  Float_t   xf[3];
  Float_t   pf[3];
  Float_t   hf[9];
  Float_t   chtr;
  Float_t   dedx2;
  Float_t   backtr;
  Float_t   cutek;
  Float_t   tlgcm2;
  Float_t   tlrad;
} Erwork_t;

//----------GCCHAN
//      COMMON/GCCHAN/LSAMVL
typedef struct {
  Int_t    lsamvl;
} Gcchan_t;

/************************************************************************
 *                                                                      *
 *      Commons for GEANE                                               *
 *                                                                      *
 ************************************************************************/

class TGeant3 : public TVirtualMC {

public:
  TGeant3();
  TGeant3(const char *title, Int_t nwgeant=0);
  virtual ~TGeant3();

  virtual void LoadAddress();
  virtual Bool_t  IsRootGeometrySupported() const {return kFALSE;}

///////////////////////////////////////////////////////////////////////
//                                                                   //
//                                                                   //
//     Here are the service routines from the geometry               //
//     which could be implemented also in other geometries           //
//                                                                   //
//                                                                   //
///////////////////////////////////////////////////////////////////////

  void  GeomIter();
  Int_t CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, 
                        Float_t &radl, Float_t &absl) const;
  Int_t NextVolUp(Text_t *name, Int_t &copy);
  Int_t CurrentVolID(Int_t &copy) const;
  Int_t CurrentVolOffID(Int_t off, Int_t &copy) const;
  const char* CurrentVolName() const;
  const char *CurrentVolOffName(Int_t off) const;
  const char* CurrentVolPath();
  Int_t VolId(const Text_t *name) const;
  Int_t MediumId(const Text_t *name) const;
  Int_t IdFromPDG(Int_t pdg) const;
  Int_t PDGFromId(Int_t pdg) const;
  const char* VolName(Int_t id) const;
  Double_t Xsec(char* reac, Double_t energy, Int_t part, Int_t mate);
  void  TrackPosition(TLorentzVector &xyz) const;
  void  TrackPosition(Double_t &x, Double_t &y, Double_t &z) const;
  void  TrackMomentum(TLorentzVector &xyz) const;
  void  TrackMomentum(Double_t &px, Double_t &py, Double_t &pz,
                      Double_t &etot) const;
  Int_t NofVolumes() const;
  Int_t NofVolDaughters(const char* volName) const;
  const char*  VolDaughterName(const char* volName, Int_t i) const;
  Int_t        VolDaughterCopyNo(const char* volName, Int_t i) const;
  Int_t    VolId2Mate(Int_t id) const;
  Double_t TrackTime() const;
  Double_t TrackCharge() const;
  Double_t TrackMass() const;
  Double_t TrackStep() const;
  Double_t TrackLength() const;
  Int_t   TrackPid() const;
  Bool_t IsNewTrack() const;
  Bool_t IsTrackInside() const;
  Bool_t IsTrackEntering() const;
  Bool_t IsTrackExiting() const;
  Bool_t IsTrackOut() const;
  Bool_t IsTrackDisappeared() const;
  Bool_t IsTrackStop() const;
  Bool_t IsTrackAlive() const;
  Int_t  NSecondaries() const;
  Int_t  CurrentEvent() const;
  TMCProcess  ProdProcess(Int_t isec) const;
  Int_t  StepProcesses(TArrayI &proc) const;
  void   GetSecondary(Int_t isec, Int_t& ipart, TLorentzVector &x,
                      TLorentzVector &p);
  Bool_t SecondariesAreOrdered() const {return kTRUE;}
  void   StopTrack();
  void   StopEvent();
  void   StopRun();
  Double_t MaxStep() const;
  void   SetMaxStep(Double_t maxstep);
  void   SetMaxNStep(Int_t maxnstp);
  Int_t  GetMaxNStep() const;
  void   ForceDecayTime(Float_t time);
  void   SetSkipNeutrinos(Bool_t flag) {fSkipNeutrinos = flag;}
  Bool_t SkipNeutrinos() {return fSkipNeutrinos;}
  Bool_t SetCut(const char* cutName, Double_t cutValue);
  Bool_t SetProcess(const char* flagName, Int_t flagValue);
  const char *GetPath();
  const char *GetNodeName();
  Bool_t DefineParticle(Int_t pdg, const char* name, 
                   TMCParticleType mcType,
                   Double_t mass, Double_t charge, Double_t lifetime);
  Bool_t DefineParticle(Int_t pdg, const char* name, 
                   TMCParticleType mcType,
                   Double_t mass, Double_t charge, Double_t lifetime,
                   const TString& /*pType*/, Double_t /*width*/, 
                   Int_t /*iSpin*/, Int_t /*iParity*/, Int_t /*iConjugation*/, 
                   Int_t /*iIsospin*/, Int_t /*iIsospinZ*/, Int_t /*gParity*/,
                   Int_t /*lepton*/, Int_t /*baryon*/,
                   Bool_t /*stable*/, Bool_t /*shortlived*/ = kFALSE,
                   const TString& /*subType*/ = "",
                   Int_t /*antiEncoding*/ = 0, Double_t /*magMoment*/ = 0.0,
                   Double_t /*excitation*/ = 0.0);
  Bool_t DefineIon(const char* name, Int_t Z, Int_t A, Int_t Q,
                   Double_t excEnergy, Double_t mass);
  virtual TString   ParticleName(Int_t pdg) const;
  virtual Double_t  ParticleMass(Int_t pdg) const;
  virtual Double_t  ParticleCharge(Int_t pdg) const;
  virtual Double_t  ParticleLifeTime(Int_t pdg) const;
  virtual TMCParticleType ParticleMCType(Int_t pdg) const;

  virtual Int_t CurrentMedium() const;
  virtual Int_t GetMedium() const;
  virtual Double_t Edep() const;
  virtual Double_t Etot() const;

  virtual void  Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
                         Double_t dens, Double_t radl, Double_t absl,
                         Float_t* buf=0, Int_t nwbuf=0);
  virtual void  Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
                         Double_t dens, Double_t radl, Double_t absl,
                         Double_t* buf, Int_t nwbuf);
  virtual void  Mixture(Int_t& kmat, const char* name, Float_t* a,Float_t* z,
                        Double_t dens, Int_t nlmat, Float_t* wmat);
  virtual void  Mixture(Int_t& kmat, const char* name, Double_t* a,Double_t* z,
                        Double_t dens, Int_t nlmat, Double_t* wmat);
  virtual void  Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
                       Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                       Double_t stemax, Double_t deemax, Double_t epsil,
                       Double_t stmin, Float_t* ubuf=0, Int_t nbuf=0);
  virtual void  Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
                       Int_t ifield, Double_t fieldm, Double_t tmaxfd,
                       Double_t stemax, Double_t deemax, Double_t epsil,
                       Double_t stmin, Double_t* ubuf, Int_t nbuf);
  virtual void Matrix(Int_t& krot, Double_t thex, Double_t phix, Double_t they,
                      Double_t phiy, Double_t thez, Double_t phiz);

  virtual void SetRootGeometry();
  virtual void SetUserParameters(Bool_t isUserParameters);

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Here are the new functions to get geometry information        //
//                        By: Bjorn S. Nilsen                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////
  // Returns the Transformation maxtrix between the volume specified by
  // the path volumePath and the Top or master volume.
  Bool_t GetTransformation(const TString &volumePath,TGeoHMatrix &mat);
  // Returns the name of the shape and its parameters for the volume
  // specified by the volumePath and the Top or master volume.
  Bool_t GetShape(const TString &volumePath,TString &shapeType,TArrayD &par);
  // Returns the material parameters for the volume specified by
  // the volume name.
  Bool_t GetMaterial(const TString &volumeName,
                     TString &name,Int_t &imat,
                     Double_t &a,Double_t &z,Double_t &den,
                     Double_t &radl,Double_t &inter,TArrayD &par);
  // Returns the medium parameters for the volume specified by the
  // volume name.
  Bool_t GetMedium(const TString &volumeName,TString &name,Int_t &imed,
                  Int_t &nmat,Int_t &isvol,Int_t &ifield,
                  Double_t &fieldm,Double_t &tmaxfd,Double_t &stemax,
                  Double_t &deemax,Double_t &epsil, Double_t &stmin,
                  TArrayD &par);

////////////////////////////////////////////////////////////////////////
//                                                                    //
//                                                                    //
//     Here are the interface functions with GEANT3.21                //
//                                                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

  // access functions to commons

   virtual Quest_t* Quest() const {return fQuest;}
   virtual Gcbank_t* Gcbank() const {return fGcbank;}
   virtual Gclink_t* Gclink() const {return fGclink;}
   virtual Gccuts_t* Gccuts() const {return fGccuts;}
   virtual Gcmore_t* Gcmore() const {return fGcmore;}
   virtual Gcmulo_t* Gcmulo() const {return fGcmulo;}
   virtual Gcmate_t* Gcmate() const {return fGcmate;}
   virtual Gctpol_t* Gctpol() const {return fGctpol;}
   virtual Gcnum_t* Gcnum() const {return fGcnum;}
   virtual Gcsets_t* Gcsets() const {return fGcsets;}
   virtual Gcopti_t* Gcopti() const {return fGcopti;}
   virtual Gctlit_t* Gctlit() const {return fGctlit;}
   virtual Gcvdma_t* Gcvdma() const {return fGcvdma;}
   virtual Gcvolu_t* Gcvolu() const {return fGcvolu;}
   virtual Gckine_t* Gckine() const {return fGckine;}
   virtual Gcflag_t* Gcflag() const {return fGcflag;}
   virtual Gctmed_t* Gctmed() const {return fGctmed;}
   virtual Gcphys_t* Gcphys() const {return fGcphys;}
   virtual Gcphlt_t* Gcphlt() const {return fGcphlt;}
   virtual Gcking_t* Gcking() const {return fGcking;}
   virtual Gckin2_t* Gckin2() const {return fGckin2;}
   virtual Gckin3_t* Gckin3() const {return fGckin3;}
   virtual Gctrak_t* Gctrak() const {return fGctrak;}
   virtual Int_t* Iq() const {return fZiq;}
   virtual Int_t* Lq() const {return fZlq;}
   virtual Float_t* Q() const {return fZq;}

  // Access to GEANE commons

   virtual Ertrio_t* Ertrio() const {return fErtrio;}
   virtual Eropts_t* Eropts() const {return fEropts;}
   virtual Eroptc_t* Eroptc() const {return fEroptc;}
   virtual Erwork_t* Erwork() const {return fErwork;}
   virtual Trcom3_t* Trcom3() const {return fTrcom3;}
   virtual Gconst_t* Gconst() const {return fGconst;}
   virtual Gconsx_t* Gconsx() const {return fGconsx;}
   virtual Gcjump_t* Gcjump() const {return fGcjump;}


      // functions from GBASE
   virtual  void  Gpcxyz();
   virtual  void  Ggclos();
   virtual  void  Gfile(const char *filename, const char *option="I");
   virtual  void  Glast();
   virtual  void  Gprint(const char *name);
   virtual  void  Grun();
   virtual  void  Gtrig();
   virtual  void  Gtrigc();
   virtual  void  Gtrigi();
   virtual  void  Gwork(Int_t nwork);
   virtual  void  Gzinit();

      // functions from GCONS
   virtual  void  Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z, 
                         Float_t &dens, Float_t &radl, Float_t &absl,
                         Float_t* ubuf, Int_t& nbuf);
   virtual  void  Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z, 
                         Double_t &dens, Double_t &radl, Double_t &absl,
                         Double_t* ubuf, Int_t& nbuf);
   virtual  void  Gfpart(Int_t ipart, char *name, Int_t &itrtyp,
                         Float_t &amass,Float_t &charge,Float_t &tlife) const;
   virtual  void  Gftmed(Int_t numed, char *name, Int_t &nmat, Int_t &isvol,
                         Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd,
                         Float_t &stemax, Float_t &deemax, Float_t &epsil,
                         Float_t &stmin, Float_t *buf=0, Int_t *nbuf=0);
   virtual  void  Gftmat(Int_t imate, Int_t ipart, char *chmeca, Int_t kdim,
                         Float_t* tkin, Float_t* value, Float_t* pcut,
                         Int_t &ixst);
   virtual  Float_t Gbrelm(Float_t z, Float_t t, Float_t cut);
   virtual  Float_t Gprelm(Float_t z, Float_t t, Float_t cut);
   virtual  void  Gmate();
   virtual  void  Gpart();
   virtual  void  Gsckov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                         Float_t *absco, Float_t *effic, Float_t *rindex);
   virtual  void  Gsdk(Int_t ipart, Float_t *bratio, Int_t *mode);
   virtual  void  Gsmate(Int_t imat, const char *name, Float_t a, Float_t z,
                         Float_t dens, Float_t radl, Float_t absl);
  virtual  void  Gfang( Float_t* p, Float_t& costh, Float_t& sinth, 
			 Float_t& cosph, Float_t& sinph, Int_t& rotate);
   virtual  void  Gsmixt(Int_t imat, const char *name, Float_t *a, Float_t *z,
                         Float_t dens, Int_t nlmat, Float_t *wmat);
   virtual  void  Gspart(Int_t ipart, const char *name, Int_t itrtyp,
                         Double_t amass, Double_t charge, Double_t tlife);
   virtual  void  Gstmed(Int_t numed,const char *name,Int_t nmat, Int_t isvol,
                         Int_t ifield, Float_t fieldm, Float_t tmaxfd,
                         Float_t stemax, Float_t deemax, Float_t epsil,
                         Float_t stmin);
   virtual  void  Gstpar(Int_t itmed, const char *param, Double_t parval);

   virtual void  SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			     Float_t *absco, Float_t *effic, Float_t *rindex);
   virtual void  SetCerenkov(Int_t itmed, Int_t npckov, Double_t *ppckov,
			     Double_t *absco, Double_t *effic, Double_t *rindex);
    
   // functions for definition of surfaces
   // and material properties for optical physics
   virtual void  DefineOpSurface(const char* name,
                         EMCOpSurfaceModel model,
                         EMCOpSurfaceType surfaceType,
                         EMCOpSurfaceFinish surfaceFinish,
                         Double_t sigmaAlpha);
   virtual void  SetBorderSurface(const char* name,
                         const char* vol1Name, int vol1CopyNo,
                         const char* vol2Name, int vol2CopyNo,
                         const char* opSurfaceName);
   virtual void  SetSkinSurface(const char* name,
                         const char* volName,
                         const char* opSurfaceName);
   virtual void  SetMaterialProperty(
                         Int_t itmed, const char* propertyName, 
                         Int_t np, Double_t* pp, Double_t* values);
   virtual void  SetMaterialProperty(
                         Int_t itmed, const char* propertyName,
                         Double_t value);
   virtual void  SetMaterialProperty(
                         const char* surfaceName, const char* propertyName, 
                         Int_t np, Double_t* pp, Double_t* values);
			 
  // functions from GKINE
   virtual  void  Gfkine(Int_t itra, Float_t *vert, Float_t *pvert,
                         Int_t &ipart, Int_t &nvert);
   virtual  void  Gfvert(Int_t nvtx,Float_t *v,Int_t &ntbeam,Int_t &nttarg,
                         Float_t &tofg);
   virtual  Int_t Gskine(Float_t *plab, Int_t ipart, Int_t nv,
                         Float_t *ubuf=0, Int_t nwbuf=0);
   virtual  Int_t Gsvert(Float_t *v, Int_t ntbeam, Int_t nttarg,
                         Float_t *ubuf=0, Int_t nwbuf=0);

      // functions from GPHYS
   virtual  void  Gphysi();

      // functions from GTRAK
   virtual  void  Gdebug();
   virtual  void  Gekbin();
   virtual  void  Gfinds();
   virtual  void  Gsking(Int_t igk);
   virtual  void  Gskpho(Int_t igk);
   virtual  void  Gsstak(Int_t iflag);
   virtual  void  Gsxyz();
   virtual  void  Gtrack();
   virtual  void  Gtreve();
   virtual  void  GtreveRoot();
   virtual  void  Grndm(Float_t *rvec, Int_t len) const;
   virtual  void  Grndmq(Int_t &is1, Int_t &is2, Int_t iseq,
                         const Text_t *chopt);

      // functions from GGEOM
   virtual  void  Gdxyz(Int_t it);
   virtual  void  Gdcxyz();

      // functions from GGEOM
   virtual  void  Gdtom(Float_t *xd, Float_t *xm, Int_t iflag);
   virtual  void  Gdtom(Double_t *xd, Double_t *xm, Int_t iflag);
   virtual  void  Glmoth(const char* iudet, Int_t iunum, Int_t &nlev,
                         Int_t *lvols, Int_t *lindx);
   virtual  void  Gmedia(Float_t *x, Int_t &numed);
   virtual  void  Gmtod(Float_t *xm, Float_t *xd, Int_t iflag);
   virtual  void  Gmtod(Double_t *xm, Double_t *xd, Int_t iflag);
   virtual  void  Gsdvn(const char *name, const char *mother,
                        Int_t ndiv, Int_t iaxis);
   virtual  void  Gsdvn2(const char *name, const char *mother,
                         Int_t ndiv, Int_t iaxis, Double_t c0i, Int_t numed);
   virtual  void  Gsdvs(const char *name, const char *mother,
                        Float_t step, Int_t iaxis, Int_t numed);
   virtual  void  Gsdvs2(const char *name, const char *mother,
                         Float_t step, Int_t iaxis, Float_t c0, Int_t numed);
   virtual  void  Gsdvt(const char *name, const char *mother,
                        Double_t step, Int_t iaxis, Int_t numed, Int_t ndvmx);
   virtual  void  Gsdvt2(const char *name, const char *mother,
                         Double_t step, Int_t iaxis,
                         Double_t c0, Int_t numed, Int_t ndvmx);
   virtual  void  Gsord(const char *name, Int_t iax);
   virtual  void  Gspos(const char *name, Int_t nr, const char *mother,
                        Double_t x, Double_t y, Double_t z, Int_t irot,
                        const char *konly="ONLY");
   virtual  void  Gsposp(const char *name, Int_t nr, const char *mother,
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np);
   virtual  void  Gsposp(const char *name, Int_t nr, const char *mother,
                         Double_t x, Double_t y, Double_t z, Int_t irot,
                         const char *konly, Double_t *upar, Int_t np);
   virtual  void  Gsrotm(Int_t nmat, Float_t theta1, Float_t phi1,
                         Float_t theta2, Float_t phi2,
                         Float_t theta3, Float_t phi3);
   virtual  void  Gprotm(Int_t nmat=0);
   virtual  Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,
                         Float_t *upar, Int_t np);
   virtual  Int_t Gsvolu(const char *name, const char *shape, Int_t nmed,
                         Double_t *upar, Int_t np);
   virtual  void Gsatt(const char *name, const char *att, Int_t val);
   virtual  void Gfpara(const char *name,Int_t number,Int_t intext,Int_t& npar,
                        Int_t& natt, Float_t* par, Float_t* att);
   virtual  void  Gckpar(Int_t ish, Int_t npar, Float_t *par);
   virtual  void  Gckmat(Int_t itmed, char *natmed);
   virtual  Int_t  Glvolu(Int_t nlev, Int_t *lnam,Int_t *lnum);
   virtual  void  Gsbool(const char* /*onlyVolName*/,
                         const char* /*manyVolName*/) {}

      // functions from GDRAW
   virtual  void  DefaultRange();
   virtual  void  InitHIGZ();
   virtual  void  Gdopen(Int_t view);
   virtual  void  Gdclose();
   virtual  void  Gdelete(Int_t view);
   virtual  void  Gdshow(Int_t view);
   virtual  void  Gdopt(const char *name,const char *value);
   virtual  void  Gdraw(const char *name,Double_t theta=30, Double_t phi=30,
                        Double_t psi=0,Double_t u0=10,Double_t v0=10,
                        Double_t ul=0.01,Double_t vl=0.01);
   virtual  void  Gdrawc(const char *name,Int_t axis=1, Float_t cut=0,
                         Float_t u0=10,Float_t v0=10,Float_t ul=0.01,
                         Float_t vl=0.01);
   virtual  void  Gdrawx(const char *name,Float_t cutthe, Float_t cutphi, 
                         Float_t cutval, Float_t theta=30, Float_t phi=30,
                         Float_t u0=10,Float_t v0=10,Float_t ul=0.01,
                         Float_t vl=0.01);
   virtual  void  Gdhead(Int_t isel, const char *name, Double_t chrsiz=0.6);
   virtual  void  Gdman(Double_t u0, Double_t v0, const char *type="MAN");
   virtual  void  Gdspec(const char *name);
   virtual  void  DrawOneSpec(const char *name);
   virtual  void  Gdtree(const char *name,Int_t levmax=15,Int_t ispec=0);
   virtual  void  GdtreeParent(const char *name,Int_t levmax=15,Int_t ispec=0);

   virtual  void  WriteEuclid(const char* filnam, const char* topvol,
			      Int_t number, Int_t nlevel);

   virtual  void  SetABAN(Int_t par=1);
   virtual  void  SetANNI(Int_t par=1);
   virtual  void  SetAUTO(Int_t par=1);
   virtual  void  SetBOMB(Float_t bomb=1);
   virtual  void  SetBREM(Int_t par=1);
   virtual  void  SetCKOV(Int_t par=1);
   virtual  void  SetClipBox(const char *name,Double_t xmin=-9999,
                             Double_t xmax=0, Double_t ymin=-9999,
                             Double_t ymax=0,Double_t zmin=-9999,
                             Double_t zmax=0);
   virtual  void  SetCOMP(Int_t par=1);
   //modified by Andrea Fontana and Alberto Rotondi - march 2007
   //added array of 5 user definable cuts (like in old Geant)
   virtual  void  SetCUTS(Float_t cutgam,Float_t cutele,Float_t cutneu,
                          Float_t cuthad,Float_t cutmuo ,Float_t bcute ,
                          Float_t bcutm ,Float_t dcute ,
                          Float_t dcutm ,Float_t ppcutm, Float_t tofmax, Float_t
			  *gcuts);
   virtual  void  InitGEANE();
   virtual void   SetClose(Int_t iclose,Float_t *pf,Float_t dstrt,
                           Float_t *w1,Float_t *w2,
			   Float_t *p1,Float_t *p2,Float_t *p3,Float_t *cl);
   virtual void   GetClose(Float_t *p1,Float_t *p2,Float_t *p3, Float_t *len);
   virtual void   SetECut(Float_t gcalpha);
   virtual  void  SetDCAY(Int_t par=1);
   virtual  void  SetDEBU(Int_t emin=1, Int_t emax=999, Int_t emod=1);
   virtual  void  SetDRAY(Int_t par=1);
   virtual  void  SetERAN(Float_t ekmin=1.e-5, Float_t ekmax=1.e4,
			  Int_t nekbin=90);
   virtual  void  SetHADR(Int_t par=1);
   virtual  void  SetKINE(Int_t kine, Float_t xk1=0, Float_t xk2=0, 
                          Float_t xk3=0, Float_t xk4=0,
                          Float_t xk5=0, Float_t xk6=0, Float_t xk7=0, 
                          Float_t xk8=0, Float_t xk9=0, Float_t xk10=0);
   virtual  void  SetLOSS(Int_t par=2);
   virtual  void  SetMULS(Int_t par=1);
   virtual  void  SetMUNU(Int_t par=1);
   virtual  void  SetOPTI(Int_t par=2);
   virtual  void  SetPAIR(Int_t par=1);
   virtual  void  SetPFIS(Int_t par=1);
   virtual  void  SetPHOT(Int_t par=1);
   virtual  void  SetRAYL(Int_t par=1);
   virtual  void  SetSTRA(Int_t par=0);
   virtual  void  SetSWIT(Int_t sw, Int_t val=1);
   virtual  void  SetTRIG(Int_t nevents=1);
   virtual  void  SetUserDecay(Int_t ipart);
   virtual  Bool_t  SetDecayMode(Int_t pdg, Float_t bratio[6], Int_t mode[6][3]);
   virtual  void  Vname(const char *name, char *vname);
   virtual  void  InitLego();

  // Routines from GEANE

   virtual void Ertrgo();
   virtual void Ertrak(const Float_t *x1, const Float_t *p1,
			const Float_t *x2, const Float_t *p2,
			Int_t ipa,  Option_t *chopt);
   virtual void Erxyzc();
   virtual void Eufill(Int_t n,Float_t *ein,Float_t *xlf);
   virtual void Eufilp(const Int_t n, Float_t *ein,
			Float_t *pli, Float_t *plf);
   virtual void Eufilv(Int_t n, Float_t *ein,
			Char_t *namv, Int_t *numv,Int_t *iovl);
   virtual void Trscsd(Float_t *pc,Float_t *rc,Float_t *pd,Float_t *rd,Float_t *h,
			Float_t ch,Int_t ierr,Float_t spu,Float_t *dj,Float_t *dk);
   virtual void Trsdsc(Float_t *pd,Float_t *rd,Float_t *pc,Float_t *rc,Float_t *h,
			Float_t *ch,Int_t *ierr,Float_t *spu,Float_t *dj,Float_t *dk);
   virtual void Trscsp(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,
			Float_t *ch,Int_t *ierr,Float_t *spx);
   virtual void Trspsc(Float_t *ps,Float_t *rs,Float_t *pc,Float_t *rc,Float_t *h,
			Float_t *ch,Int_t *ierr,Float_t *spx);


  // Control Methods

  virtual void FinishGeometry();
  virtual void BuildPhysics();
  virtual void Init();
  virtual void ProcessEvent();
  virtual Bool_t ProcessRun(Int_t nevent);
  virtual void AddParticlesToPdgDataBase() const;

  //
  virtual void SetColors();

  void  SetTrack(Int_t done, Int_t parent, Int_t pdg,
  	          Float_t *pmom, Float_t *vpos, Float_t *polar,
                  Float_t tof, TMCProcess mech, Int_t &ntr,
                  Float_t weight, Int_t is);

  Ertrio_t  *fErtrio;          //! ERTRIO common structure
  Ertrio1_t *fErtrio1;         //! ERTRIO1 common structure
  Eropts_t  *fEropts;          //! EROPTS common structure
  Eroptc_t  *fEroptc;          //! EROPTC common structure
  Erwork_t  *fErwork;          //! ERWORK common structure
  Trcom3_t  *fTrcom3;          //! TRCOM3 common structure


private:
  Int_t ConvertVolumePathString(const TString &volumeName,Int_t **lnam,
                                Int_t **lnum);



protected:
  Int_t fNextVol;    // Iterator for GeomIter
  char  fPath[512];  // Current path of G3
//--------------Declarations for ZEBRA---------------------
  Int_t *fZiq;                //! Good Old IQ of Zebra
  Int_t *fZlq;                //! Good Old LQ of Zebra
  Float_t *fZq;               //! Good Old Q of Zebra

  Quest_t  *fQuest;           //! QUEST common structure
  Gcbank_t *fGcbank;          //! GCBANK common structure
  Gclink_t *fGclink;          //! GCLINK common structure
  Gccuts_t *fGccuts;          //! GCCUTS common structure
  Gcmore_t *fGcmore;          //! GCMORE common structure
  Gcmulo_t *fGcmulo;          //! GCMULO common structure
  Gcmate_t *fGcmate;          //! GCMATE common structure
  Gctpol_t *fGctpol;          //! GCTPOL common structure
  Gcnum_t  *fGcnum;           //! GCNUM common structure
  Gcsets_t *fGcsets;          //! GCSETS common structure
  Gcopti_t *fGcopti;          //! GCOPTI common structure
  Gctlit_t *fGctlit;          //! GCTLIT common structure
  Gcvdma_t *fGcvdma;          //! GCVDMA common structure
  Gcvolu_t *fGcvolu;          //! GCVOLU common structure
  Gckine_t *fGckine;          //! GCKINE common structure
  Gcflag_t *fGcflag;          //! GCFLAG common structure
  Gctmed_t *fGctmed;          //! GCTMED common structure
  Gcphys_t *fGcphys;          //! GCPHYS common structure
  Gcphlt_t *fGcphlt;          //! GCPHLT common structure
  Gcking_t *fGcking;          //! GCKING common structure
  Gckin2_t *fGckin2;          //! GCKIN2 common structure
  Gckin3_t *fGckin3;          //! GCKIN3 common structure
  Gctrak_t *fGctrak;          //! GCTRAK common structure
  Gcchan_t *fGcchan;          //! GCCHAN common structure


  // commons for GEANE
  Gconst_t *fGconst;          //! GCONST common structure
  Gconsx_t *fGconsx;          //! GCONSX common structure
  Gcjump_t *fGcjump;          //! GCJUMP common structure


 

  //Put here all volume names

  char (*fVolNames)[5];           //! Names of geant volumes as C++ chars
  TObjArray fMedNames;            //! Names of geant medias as TObjString

  Int_t fNG3Particles;            // Number of G3 particles
  Int_t fNPDGCodes;               // Number of PDG codes known by G3

  TArrayI          fPDGCode;// Translation table of PDG codes
  TGeoMCGeometry*  fMCGeo;  // Implementation of TVirtualMCGeometry for TGeo
  Bool_t           fImportRootGeometry; // Option to import geometry from TGeo
                                        // (materials and medias are filled
                                        // in FinishGeometry()
  Bool_t           fStopRun;            // The flag for stopping run by a user
  Bool_t           fSkipNeutrinos;      // The flag for skipping neutrinos from decays

  TMCProcess G3toVMC(Int_t iproc) const;

  void   DefineParticles();
  Int_t  TransportMethod(TMCParticleType particleType) const;
  TString  ParticleClass(TMCParticleType particleType) const;
  TMCParticleType ParticleType(Int_t itrtyp) const;

  enum {kTRIG = BIT(14),
        kSWIT = BIT(15),
        kDEBU = BIT(16),
        kAUTO = BIT(17),
        kABAN = BIT(18),
        kOPTI = BIT(19),
        kERAN = BIT(20)
  };
  TGeant3(const TGeant3 &) : TVirtualMC() {}
  TGeant3 & operator=(const TGeant3&) {return *this;}

  // array conversion
  Float_t* CreateFloatArray(Float_t*  array, Int_t size) const;
  Float_t* CreateFloatArray(Double_t* array, Int_t size) const;
  Int_t    NextKmat() const;

  // functions for building geometry with different interface
  // for double and single precision
  void  G3Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
	            Double_t dens, Double_t radl, Double_t absl,
	            Float_t* buf=0, Int_t nwbuf=0);
  void  G3Mixture(Int_t& kmat, const char* name, Float_t* a,Float_t* z,
	            Double_t dens, Int_t nlmat, Float_t* wmat);
  void  G3Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
	            Int_t ifield, Double_t fieldm, Double_t tmaxfd,
		    Double_t stemax, Double_t deemax, Double_t epsil,
	            Double_t stmin, Float_t* ubuf=0, Int_t nbuf=0);
  Int_t G3Gsvolu(const char *name, const char *shape, Int_t nmed,
                    Float_t *upar, Int_t np);
  void  G3Gsposp(const char *name, Int_t nr, const char *mother,
                 Double_t x, Double_t y, Double_t z, Int_t irot, 
                 const char *konly, Float_t *upar, Int_t np);
                 
  // particles definition
  Int_t GetIonPdg(Int_t z, Int_t a, Int_t i = 0) const;                
  Int_t GetSpecialPdg(Int_t number) const;                

  ClassDef(TGeant3,1)  //C++ interface to Geant basic routines
};

#endif //ROOT_TGeant3
