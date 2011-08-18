#include <iostream>
#include "TGeant3.h"
#include "TCallf77.h" 
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TParticle.h"

//#define COLLECT_TRACKS
#if defined(COLLECT_TRACKS)
#include "TGeoManager.h"
#include "TVirtualGeoTrack.h"
#endif
   
#ifndef WIN32
#  define gudigi gudigi_
#  define guhadr guhadr_
#  define guout  guout_
#  define guphad guphad_
#  define gudcay gudcay_
#  define guiget guiget_
#  define guinme guinme_
#  define guinti guinti_
#  define gunear gunear_
#  define guskip guskip_
#  define guview guview_
#  define gupara gupara_
#  define gudtim gudtim_
#  define guplsh guplsh_
#  define gutrev gutrev_
#  define gutrak gutrak_
#  define guswim guswim_
#  define gufld  gufld_
#  define gustep gustep_
#  define gukine gukine_

#  define eustep eustep_

#  define calsig calsig_
#  define gcalor gcalor_

#  define gheish gheish_
#  define flufin flufin_
#  define gfmfin gfmfin_
#  define gpghei gpghei_
#  define fldist fldist_
#  define gfmdis gfmdis_
#  define g3helx3 g3helx3_
#  define g3helix g3helix_
#  define g3rkuta g3rkuta_
#  define g3track g3track_
#  define gtreveroot gtreveroot_
#  define g3last  g3last_
#  define g3invol g3invol_
#  define g3tmedi g3tmedi_
#  define g3media g3media_
#  define g3tmany g3tmany_
#  define g3tnext g3tnext_
#  define ginvol ginvol_
#  define gtmedi gtmedi_
#  define gtmany gtmany_
#  define gtonly gtonly_
#  define gmedia gmedia_
#  define glvolu glvolu_
#  define gtnext gtnext_

#else
#  define gudigi GUDIGI
#  define guhadr GUHADR
#  define guout  GUOUT
#  define guphad GUPHAD
#  define gudcay GUDCAY
#  define guiget GUIGET
#  define guinme GUINME
#  define guinti GUINTI
#  define gunear GUNEAR
#  define guskip GUSKIP
#  define guview GUVIEW
#  define gupara GUPARA
#  define gudtim GUDTIM
#  define guplsh GUPLSH
#  define gutrev GUTREV
#  define gutrak GUTRAK
#  define guswim GUSWIM
#  define gufld  GUFLD
#  define gustep GUSTEP
#  define gukine GUKINE

#  define calsig CALSIG
#  define gcalor GCALOR

#  define gheish GHEISH
#  define flufin FLUFIN
#  define gfmfin GFMFIN
#  define gpghei GPGHEI
#  define fldist FLDIST
#  define gfmdis GFMDIS
#  define g3helx3 G3HELX3
#  define g3helix G3HELIX
#  define g3rkuta G3RKUTA
#  define gtrack GTRACK
#  define gtreveroot GTREVEROOT
#  define glast  GLAST
#  define ginvol GINVOL
#  define gtmedi GTMEDI
#  define gtmany GTMANY
#  define gmedia GMEDIA
#  define glvolu GLVOLU
#  define gtnext GTNEXT 

#endif

extern TGeant3* geant3;
#if defined(COLLECT_TRACKS)
extern TGeoManager *gGeoManager;
#endif


extern "C" type_of_call void calsig();
extern "C" type_of_call void gcalor();

extern "C" type_of_call void gheish();
extern "C" type_of_call void flufin();
extern "C" type_of_call void gfmfin();
extern "C" type_of_call void gpghei();
extern "C" type_of_call void fldist();
extern "C" type_of_call void gfmdis();
extern "C" type_of_call void g3helx3(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void g3helix(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void g3rkuta(Float_t&, Float_t&, Float_t*, Float_t*);
extern "C" type_of_call void g3track();
extern "C" type_of_call void gtreveroot();
extern "C" type_of_call void g3last();
extern "C" type_of_call void g3invol(Float_t*, Int_t&);
extern "C" type_of_call void g3tmedi(Float_t*, Int_t&);
extern "C" type_of_call void g3tmany(Int_t&);
extern "C" type_of_call void g3media(Float_t*, Int_t&, Int_t&);
extern "C" type_of_call void g3tnext();
extern "C" type_of_call void ginvol(Float_t*, Int_t&);
extern "C" type_of_call void gtmedi(Float_t*, Int_t&);
extern "C" type_of_call void gtmany(Int_t&);
extern "C" type_of_call void gtonly(Int_t&);
extern "C" type_of_call void gmedia(Float_t*, Int_t&, Int_t&);
extern "C" type_of_call void glvolu(Int_t &nlev, Int_t *lnam,Int_t *lnum, Int_t &ier);
extern "C" type_of_call void gtnext();
extern "C" type_of_call {

//______________________________________________________________________
void gudigi() 
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to digitize one event                       *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************

}


//______________________________________________________________________
void guhadr()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to generate one hadronic interaction        *
//    *                                                                *
//    *    ==>Called by : GTHADR,GTNEUT                                *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      Int_t ihadr = geant3->Gcphys()->ihadr;

      //    printf(" iHadr 1 :%i \n ",ihadr);
      if (ihadr<4)      { gheish();}
      else if (ihadr==4){ flufin();}
      else if (ihadr==5){ gcalor();}
      else              { gfmfin();}
      
}

//______________________________________________________________________
void guout()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called at the end of each event             *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
  //printf("count_gmedia= %8d\n",count_gmedia);
  //printf("count_gtmedi= %8d\n",count_gtmedi);
  //printf("count_ginvol= %8d\n",count_ginvol);
  //printf("count_gtnext= %8d\n",count_gtnext);
}

//______________________________________________________________________
void guphad()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to compute Hadron. inter. probabilities     *
//    *                                                                *
//    *    ==>Called by : GTHADR,GTNEUT                                *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      Int_t ihadr = geant3->Gcphys()->ihadr;
      // printf(" iHadr 2 :%i \n ",ihadr);
      if (ihadr<4){ gpghei();}
      else if ( ihadr==4 ){ fldist();}
      else if ( ihadr==5 ){ calsig();}
      else    { gfmdis();}

}

//______________________________________________________________________
void gudcay()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to decay particles                          *
//    *                                                                *
//    *    ==>Called by : G3DECAY                                      *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
    
    // set decay table
    if (!gMC->GetDecayer()) return;
    gMC->GetDecayer()->ForceDecay();

// Initialize 4-momentum vector    
    Int_t ipart = geant3->Gckine()->ipart;
    static TLorentzVector p;

    Float_t pmom = geant3->Gctrak()->vect[6];

    p[0] = pmom * (geant3->Gctrak()->vect[3]);
    p[1] = pmom * (geant3->Gctrak()->vect[4]);
    p[2] = pmom * (geant3->Gctrak()->vect[5]);
    p[3] = geant3->Gctrak()->getot;
    
    
// Convert from geant to lund particle code
    Int_t iplund=gMC->PDGFromId(ipart);
    
// Particle list
    static TClonesArray *particles;
    if(!particles) particles=new TClonesArray("TParticle",1000);
// Decay
    gMC->GetDecayer()->Decay(iplund, &p);
    
// Fetch Particles
    Int_t np = geant3->GetDecayer()->ImportParticles(particles);
    if (np <=1) return;

    TParticle *  iparticle = (TParticle *) particles->At(0);
    Int_t ipF = 0, ipL = 0 ;
    Int_t i,j;

// Array to flag deselected particles
    Int_t*  pFlag = new Int_t[np];
    for (i=0; i<np; i++) pFlag[i]=0;
// Particle loop
    for (i=1; i < np; i++) 
    {
	iparticle = (TParticle *) particles->At(i);
	ipF = iparticle->GetFirstDaughter();
	ipL = iparticle->GetLastDaughter();	
	Int_t kf = iparticle->GetPdgCode();
	Int_t ks = iparticle->GetStatusCode();
//
// Deselect daughters of deselected particles
// and jump skip the current particle
	if (pFlag[i] == 1) {
	    if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
	    continue;
	} // deselected ??
// Particles with long life-time are put on the stack for further tracking
// Decay products are deselected
//	
	if (ks != 1) { 
	    Double_t lifeTime = gMC->GetDecayer()->GetLifetime(kf);
	    if (lifeTime > (Double_t) 1.e-15) {
		if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
	    } else{
		continue;
	    }
	} // ks==1 ?
// Skip neutrinos
        if (geant3->SkipNeutrinos()) {
           if (kf==12 || kf ==-12) continue;
           if (kf==14 || kf ==-14) continue;
           if (kf==16 || kf ==-16) continue;
        }
	
	Int_t index=geant3->Gcking()->ngkine;
// Put particle on geant stack
// momentum vector
	
	(geant3->Gcking()->gkin[index][0]) = iparticle->Px();
	(geant3->Gcking()->gkin[index][1]) = iparticle->Py();
	(geant3->Gcking()->gkin[index][2]) = iparticle->Pz();
	(geant3->Gcking()->gkin[index][3]) = iparticle->Energy();
	Int_t ilu = gMC->IdFromPDG(kf);

// particle type	
	(geant3->Gcking()->gkin[index][4]) = Float_t(ilu);
// position
	(geant3->Gckin3()->gpos[index][0]) = geant3->Gctrak()->vect[0];
	(geant3->Gckin3()->gpos[index][1]) = geant3->Gctrak()->vect[1];
	(geant3->Gckin3()->gpos[index][2]) = geant3->Gctrak()->vect[2];
// time of flight offset (mm)
	(geant3->Gcking()->tofd[index])    = 0.;
// increase stack counter
	(geant3->Gcking()->ngkine)=index+1;
    }
    delete[] pFlag;
}

//______________________________________________________________________
void guiget(Int_t&, Int_t&, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine for interactive control of GEANT            *
//    *                                                                *
//    *    ==>Called by : <GXINT>, GINCOM                              *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void guinme(Float_t*, Int_t&, Float_t*, Int_t& IYES)
{
//
//    **********************************************
//    *                                            *
//    *    USER ROUTINE TO PROVIDE GINME FUNCTION  *
//    *    FOR ALL USER SHAPES IDENTIFIED BY THE   *
//    *    SHAPE NUMBER SH. POINT IS GIVEN IN X    *
//    *    THE PARAMETERS ARE GIVEN IN P. IYES IS  *
//    *    RETURNED 1 IF POINT IS IN, 0 IF POINT   *
//    *    IS OUT AND LESS THAN ZERO IF SHAPE      *
//    *    NUMBER IS NOT SUPPORTED.                *
//    *                                            *
//    *    ==>Called by : GINME                    *
//    *                                            *
//    **********************************************
//
      IYES=-1;
}

//______________________________________________________________________
void guinti()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine for interactive version                     *
//    *                                                                *
//    *    ==>Called by : <GXINT>,  GINTRI                             *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void gunear(Int_t&, Int_t&, Float_t*, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *    User search                                                 *
//    *       ISEARC to identify the given volume                      *
//    *       ICALL  to identify the calling routine                   *
//    *              1 GMEDIA like                                     *
//    *              2 GNEXT like                                      *
//    *       X      coordinates (+direction for ICALL=2)              *
//    *       JNEAR  address of default list of neighbours             *
//    *              (list to be overwriten by user)                   *
//    *                                                                *
//    *    Called by : GFTRAC, GINVOL, GTMEDI, GTNEXT, GNEXT, GMEDIA   *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void guskip(Int_t& ISKIP)
{
//
//    ******************************************************************
//    *                                                                *
//    *   User routine to skip unwanted tracks                         *
//    *                                                                *
//    *   Called by : GSSTAK                                           *
//    *   Author    : F.Bruyant                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      ISKIP = 0;
}

//______________________________________________________________________
void guswim(Float_t& CHARGE, Float_t& STEP, Float_t* VECT, Float_t* VOUT)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to control tracking of one track            *
//    *       in a magnetic field                                      *
//    *                                                                *
//    *    ==>Called by : GTELEC,GTHADR,GTMUON                         *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
  Int_t   ifield = geant3->Gctmed()->ifield;
  Float_t fieldm = geant3->Gctmed()->fieldm;

  if (ifield==3) {
    Float_t fldcharge = fieldm*CHARGE;
    g3helx3(fldcharge,STEP,VECT,VOUT);
  }
  else if (ifield==2) g3helix(CHARGE,STEP,VECT,VOUT);
  else                g3rkuta(CHARGE,STEP,VECT,VOUT);
}

//______________________________________________________________________
void guview(Int_t&, Int_t&, DEFCHARD, Int_t& DEFCHARL)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine for interactive version                     *
//    *                                                                *
//    *    ==>Called by : <GXINT>, GINC1                               *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
void gupara()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called every time a particle falls below    *
//    *       parametrization threshold. This routine should create    *
//    *       the parametrization stack, and, when this is full,       *
//    *       parametrize the shower and track the geantinos.          *
//    *                                                                *
//    *    ==>Called by : GTRACK                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
}

//______________________________________________________________________
Float_t gudtim(Float_t&, Float_t&, Int_t&, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *       User function called by GCDRIF to return drift time      *
//    *                                                                *
//    *    ==>Called by : GCDRIF                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
      return 0;
}


//______________________________________________________________________
Float_t guplsh(Int_t&, Int_t&)
{
//
//    ******************************************************************
//    *                                                                *
//    *                                                                *
//    *    ==>Called by : GLISUR                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
//
//*** By default this defines perfect smoothness
      return 1;
}

//______________________________________________________________________
void gutrak()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to control tracking of one track            *
//    *                                                                *
//    *    ==>Called by : GTREVE                                       *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
     TVirtualMCApplication::Instance()->PreTrack();

     g3track();

     TVirtualMCApplication::Instance()->PostTrack();
}

//______________________________________________________________________
void gutrev()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine to control tracking of one event            *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//
  gtreveroot();
}

//______________________________________________________________________
#ifdef SINGLEFIELD
void gufld(Float_t *x, Float_t *b)
{
  Double_t xdouble[3];
  Double_t bdouble[3];
  for (Int_t i=0; i<3; i++) xdouble[i] = x[i]; 
#else
void gufld(Double_t *xdouble, Double_t *bdouble)
{
#endif
  if ( gMC->GetMagField() ) {
    gMC->GetMagField()->Field(xdouble,bdouble);
  }
  else {  
    static Bool_t warn = true;   
    if (warn) { 
      Warning("gufld", "Using deprecated function TVirtualMCApplication::Field().");
      Warning("gufld", "New TVirtualMagField interface should be used instead.");
      warn = false;
    }        

    TVirtualMCApplication::Instance()->Field(xdouble,bdouble);
  }  

#ifdef SINGLEFIELD
  for (Int_t j=0; j<3; j++) b[j] = bdouble[j]; 
#endif
}
//______________________________________________________________________
void eustep(){
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called at the end of each tracking step     *
//    *       when using GEANE                                         *
//    ******************************************************************
  
   Int_t cflag;
   static  Int_t icc=0;
   static  Int_t icont=0;
   Float_t  dist2;
   static Float_t prdist2;
   Float_t  d2x,d2y,d2z,amodd;
// clear icc when tracking starts
   if(geant3->Gctrak()->sleng == 0) icc = 0;
   cflag=geant3->Gcmore()->iclose;
   if (geant3->Gcflag()->idebug * geant3->Gcflag()->iswit[2] != 0)geant3->Erxyzc();
   if(cflag==1){
// distance between the track point and the point        
      if(icc==0) prdist2=geant3->Gconst()->big;
         dist2 = (geant3->Gctrak()->vect[0] - geant3->Gcmore()->pfinal[0])*
              (geant3->Gctrak()->vect[0] - geant3->Gcmore()->pfinal[0])+
              (geant3->Gctrak()->vect[1] - geant3->Gcmore()->pfinal[1])*
              (geant3->Gctrak()->vect[1] - geant3->Gcmore()->pfinal[1])+
              (geant3->Gctrak()->vect[2] - geant3->Gcmore()->pfinal[2])*
              (geant3->Gctrak()->vect[2] - geant3->Gcmore()->pfinal[2]);

      if((TMath::Sqrt(dist2)-TMath::Sqrt(prdist2)) < 1.e-3){
         prdist2 = dist2;
         icc = 1;
         icont = 1;
         geant3->Gcmore()->cleng[0] = geant3->Gcmore()->cleng[1];
         geant3->Gcmore()->cleng[1] = geant3->Gcmore()->cleng[2];
         geant3->Gcmore()->cleng[2] = geant3->Gctrak()->sleng;
         for(Int_t i=0; i<3; i++)geant3->Gcmore()->p1[i] = geant3->Gcmore()->p2[i];  //call ucopy(p2,p1,3)
         for(Int_t i=0; i<3; i++)geant3->Gcmore()->p2[i] = geant3->Gcmore()->p3[i];  //call ucopy(p3,p2,3)
         for(Int_t i=0; i<3; i++)geant3->Gcmore()->p3[i] = geant3->Gctrak()->vect[i]; //call ucopy(vect,p3,3)
      }else{ // store the first point of increasing distance
         if(icont == 1) {
            geant3->Gcmore()->cleng[0] = geant3->Gcmore()->cleng[1];
            geant3->Gcmore()->cleng[1] = geant3->Gcmore()->cleng[2];   
            geant3->Gcmore()->cleng[2] = geant3->Gctrak()->sleng;
            for(Int_t i=0; i<3; i++)geant3->Gcmore()->p1[i] = geant3->Gcmore()->p2[i];  //call ucopy(p2,p1,3)
            for(Int_t i=0; i<3; i++)geant3->Gcmore()->p2[i] = geant3->Gcmore()->p3[i];  //call ucopy(p3,p2,3)
            for(Int_t i=0; i<3; i++)geant3->Gcmore()->p3[i] = geant3->Gctrak()->vect[i]; //call ucopy(vect,p3,3)
            icont = 0;
         }
      } 
   }else if(cflag==2) {
   //   printf("geant3->Gconst()->big = %F" ,geant3->Gconst());
      if(icc == 0) prdist2=geant3->Gconst()->big;
      d2x = (geant3->Gcmore()->wire2[1]-geant3->Gcmore()->wire1[1])*
            (geant3->Gcmore()->wire1[2]-geant3->Gctrak()->vect[2])-
            (geant3->Gcmore()->wire2[2]-geant3->Gcmore()->wire1[2])*
            (geant3->Gcmore()->wire1[1]-geant3->Gctrak()->vect[1]); 

      d2y = (geant3->Gcmore()->wire2[2]-geant3->Gcmore()->wire1[2])*
            (geant3->Gcmore()->wire1[0]-geant3->Gctrak()->vect[0])- 
            (geant3->Gcmore()->wire2[0]-geant3->Gcmore()->wire1[0])*
            (geant3->Gcmore()->wire1[2]-geant3->Gctrak()->vect[2]);

      d2z = (geant3->Gcmore()->wire2[0]-geant3->Gcmore()->wire1[0])*
            (geant3->Gcmore()->wire1[1]-geant3->Gctrak()->vect[1])-
            (geant3->Gcmore()->wire2[1]-geant3->Gcmore()->wire1[1])*
            (geant3->Gcmore()->wire1[0]-geant3->Gctrak()->vect[0]);

      amodd =(geant3->Gcmore()->wire2[0]-geant3->Gcmore()->wire1[0])*
             (geant3->Gcmore()->wire2[0]-geant3->Gcmore()->wire1[0])+ 
             (geant3->Gcmore()->wire2[1]-geant3->Gcmore()->wire1[1])*
             (geant3->Gcmore()->wire2[1]-geant3->Gcmore()->wire1[1])+
             (geant3->Gcmore()->wire2[2]-geant3->Gcmore()->wire1[2])*
             (geant3->Gcmore()->wire2[2]-geant3->Gcmore()->wire1[2]);   

      dist2 = (d2x*d2x + d2y*d2y + d2z*d2z)/amodd;

//   distance between the track point and the wire   
      if((TMath::Sqrt(dist2)-TMath::Sqrt(prdist2)) < 1.e-3) {
         prdist2 = dist2;
         icc=1;
         icont = 1;
         geant3->Gcmore()->cleng[0] = geant3->Gcmore()->cleng[1];
         geant3->Gcmore()->cleng[1] = geant3->Gcmore()->cleng[2];
         geant3->Gcmore()->cleng[2] = geant3->Gctrak()->sleng;
         for(Int_t i=0; i<3; i++)geant3->Gcmore()->p1[i] = geant3->Gcmore()->p2[i];  //call ucopy(p2,p1,3)
         for(Int_t i=0; i<3; i++)geant3->Gcmore()->p2[i] = geant3->Gcmore()->p3[i];  //call ucopy(p3,p2,3)
         for(Int_t i=0; i<3; i++)geant3->Gcmore()->p3[i] = geant3->Gctrak()->vect[i]; //call ucopy(vect,p3,3)
      }else{ //     store the first point of increasing distance
         if(icont == 1){
            geant3->Gcmore()->cleng[0] = geant3->Gcmore()->cleng[1];
            geant3->Gcmore()->cleng[1] = geant3->Gcmore()->cleng[2];
            geant3->Gcmore()->cleng[2] = geant3->Gctrak()->sleng;
            for(Int_t i=0; i<3; i++)geant3->Gcmore()->p1[i] = geant3->Gcmore()->p2[i];  //call ucopy(p2,p1,3)
            for(Int_t i=0; i<3; i++)geant3->Gcmore()->p2[i] = geant3->Gcmore()->p3[i];  //call ucopy(p3,p2,3)
            for(Int_t i=0; i<3; i++)geant3->Gcmore()->p3[i] = geant3->Gctrak()->vect[i]; //call ucopy(vect,p3,3)
            icont = 0;
         }
      }
   }
   
    // --- Particle leaving the setup ?
   if (!gMC->IsTrackOut()) {
      TVirtualMCApplication *app = TVirtualMCApplication::Instance();
      app->GeaneStepping();
   }
   
}
//______________________________________________________________________
void gustep()
{
//
//    ******************************************************************
//    *                                                                *
//    *       User routine called at the end of each tracking step     *
//    *       INWVOL is different from 0 when the track has reached    *
//    *              a volume boundary                                 *
//    *       ISTOP is different from 0 if the track has stopped       *
//    *                                                                *
//    *    ==>Called by : GTRACK                                       *
//    *                                                                *
//    ******************************************************************
//
  Int_t ipp, jk, nt;
  Float_t polar[3]={0,0,0};
  Float_t mom[3];
  static TMCProcess pProc;
  
  TVirtualMCApplication *app = TVirtualMCApplication::Instance();
  TVirtualMCStack* stack = gMC->GetStack();
  //     Stop particle if outside user defined tracking region 
  Double_t x, y, z, rmax;
  geant3->TrackPosition(x,y,z);

#if defined(COLLECT_TRACKS)
  if (gMC->IsRootGeometrySupported()) {
  Int_t nstep = geant3->Gctrak()->nstep;
  Int_t cpdg = gMC->PDGFromId(geant3->Gckine()->ipart);
  Bool_t isnew = kFALSE; // gMC->IsNewTrack() returns true just for new used indices
  if (nstep==0) isnew = kTRUE;
  Int_t cid = stack->GetCurrentTrackNumber();
  Int_t mid = stack->GetCurrentParentTrackNumber();
  Double_t tofg = geant3->Gctrak()->tofg;
  //printf("id=%i mid=%i pdg=%i nstep=%i (%f,%f,%f)\n",cid,mid,cpdg,nstep,x,y,z);

  TVirtualGeoTrack *parent = 0;
  if (mid>=0) {
     parent = gGeoManager->GetTrackOfId(mid);
     if (!parent) printf("Error: no parent track with id=%i\n",mid);
  }   
  TVirtualGeoTrack *track;
  if (isnew) {
     if (parent) {
        //printf("Adding daughter %i of %i\n",cid,mid);
        track = parent->AddDaughter(cid, cpdg);
        gGeoManager->SetCurrentTrack(track);
     } else {     
        Int_t itrack = gGeoManager->AddTrack(cid, cpdg); 
        //printf("Added new primary %i\n",cid);
        gGeoManager->SetCurrentTrack(itrack);
        track = gGeoManager->GetCurrentTrack();
     }
     TDatabasePDG *pdgdb = TDatabasePDG::Instance();
     if (pdgdb) {
        TParticlePDG *part = pdgdb->GetParticle(cpdg);
        if (part) {
           track->SetName(part->GetName());
           track->SetParticle(part);
        }   
     }   
  } else {
     track = gGeoManager->GetCurrentTrack();
  } 
  Double_t xo,yo,zo,to;
  Bool_t skippoint = kFALSE;
  if (track->HasPoints()) {
     track->GetLastPoint(xo,yo,zo,to);
     Double_t rdist = TMath::Sqrt((xo-x)*(xo-x)+(yo-y)*(yo-y)+(zo-z)*(zo-z));
     if (rdist<0.01) skippoint=kTRUE;
  }   
  if (!skippoint) track->AddPoint(x,y,z,tofg);
  }
#endif  

  rmax = app->TrackingRmax();
  if (x*x+y*y > rmax*rmax ||
      TMath::Abs(z) > app->TrackingZmax()) {
	gMC->StopTrack();
  }

  // --- Add new created particles 
  if (gMC->NSecondaries() > 0) {
    pProc=gMC->ProdProcess(0);
    for (jk = 0; jk < geant3->Gcking()->ngkine; ++jk) {
      ipp = Int_t (geant3->Gcking()->gkin[jk][4]+0.5);
      // --- Skip neutrinos! 
      if (ipp != 4 || !(geant3->SkipNeutrinos())) {
        geant3->SetTrack(1,stack->GetCurrentTrackNumber(),gMC->PDGFromId(ipp), geant3->Gcking()->gkin[jk], 
			 geant3->Gckin3()->gpos[jk], polar,geant3->Gctrak()->tofg, pProc, nt, 1., 0);
      }
    }
  }
  // Cherenkov photons here
  if ( geant3->Gckin2()->ngphot ) {
    for (jk = 0; jk < geant3->Gckin2()->ngphot; ++jk) {
      mom[0]=geant3->Gckin2()->xphot[jk][3]*geant3->Gckin2()->xphot[jk][6];
      mom[1]=geant3->Gckin2()->xphot[jk][4]*geant3->Gckin2()->xphot[jk][6];
      mom[2]=geant3->Gckin2()->xphot[jk][5]*geant3->Gckin2()->xphot[jk][6];
      geant3->SetTrack(1, stack->GetCurrentTrackNumber(), gMC->PDGFromId(50),
		       mom,                             //momentum
		       geant3->Gckin2()->xphot[jk],     //position
		       &geant3->Gckin2()->xphot[jk][7], //polarisation
		       geant3->Gckin2()->xphot[jk][10], //time of flight
		       kPCerenkov, nt, 1., 0);
      }
  }
  // --- Particle leaving the setup ?
  if (!gMC->IsTrackOut()) app->Stepping();

  // --- Standard GEANT debug routine 
  //g3pcxyz();
  if(geant3->Gcflag()->idebug) geant3->Gdebug();
}

//______________________________________________________________________
void gukine ()
{
//
//    ******************************************************************
//    *                                                                *
//    *       Read or Generates Kinematics for primary tracks          *
//    *                                                                *
//    *    ==>Called by : GTRIG                                        *
//    *                                                                *
//    ******************************************************************
//
//
//    ------------------------------------------------------------------
//

  TVirtualMCApplication::Instance()->GeneratePrimaries();
}
}
// end of extern "C"
