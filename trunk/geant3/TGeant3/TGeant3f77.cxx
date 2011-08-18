/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Interface Class to the Geant3.21 MonteCarlo                              //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TGeant3f77.h" 

#include "TCallf77.h" 


#ifndef WIN32 
# define gzebra  gzebra_ 
# define grfile  grfile_ 
# define gpcxyz  gpcxyz_ 
# define ggclos  ggclos_ 
# define glast   glast_ 
# define ginit   ginit_ 
# define g3cinit  g3cinit_ 
# define gzinit  gzinit_ 
# define grun    grun_ 
# define gtrig   gtrig_ 
# define gtrigc  gtrigc_ 
# define gtrigi  gtrigi_ 
# define gfmate  gfmate_ 
# define gfpart  gfpart_ 
# define gftmed  gftmed_ 
# define gftmat  gftmat_ 
# define gpart   gpart_ 
# define gmate   gmate_ 
# define gsdk    gsdk_ 
# define gsmate  gsmate_ 
# define gsmixt  gsmixt_ 
# define gspart  gspart_ 
# define gstmed  gstmed_ 
# define gsckov  gsckov_
# define gstpar  gstpar_ 
# define gfkine  gfkine_ 
# define gfvert  gfvert_ 
# define gskine  gskine_ 
# define gsvert  gsvert_ 
# define gpvolu  gpvolu_ 
# define gprotm  gprotm_ 
# define gptmed  gptmed_ 
# define gpmate  gpmate_ 
# define gppart  gppart_ 
# define gpsets  gpsets_ 
# define gpvert  gpvert_ 
# define gpkine  gpkine_ 
# define gpjxyz  gpjxyz_ 
# define gphits  gphits_ 
# define g3pvolu  g3pvolu_ 
# define g3protm  g3protm_ 
# define g3ptmed  g3ptmed_ 
# define g3pmate  g3pmate_ 
# define g3ppart  g3ppart_ 
# define g3psets  g3psets_ 
# define g3pvert  g3pvert_ 
# define g3pkine  g3pkine_ 
# define g3pjxyz  g3pjxyz_ 
# define g3phits  g3phits_ 
# define g3part  g3part_ 
# define g3mate  g3mate_ 
# define gscank  gscank_ 
# define gscanu  gscanu_ 
# define g3scank  g3scank_ 
# define g3scanu  g3scanu_ 
# define g3bhsta  g3bhsta_ 
# define gbhsta  gbhsta_ 
# define gphysi  gphysi_ 
# define gdebug  gdebug_ 
# define gekbin  gekbin_ 
# define gfinds  gfinds_ 
# define gsking  gsking_ 
# define gskpho  gskpho_ 
# define gsstak  gsstak_ 
# define gsxyz   gsxyz_ 
# define gtrack  gtrack_ 
# define gtreve  gtreve_ 
# define gtreveroot  gtreveroot_ 
# define grndm   grndm_ 
# define grndmq  grndmq_ 
# define gdtom   gdtom_ 
# define glmoth  glmoth_ 
# define gmtod   gmtod_ 
# define gsdvn   gsdvn_ 
# define gsdvn2  gsdvn2_ 
# define gsdvs   gsdvs_ 
# define gsdvs2  gsdvs2_ 
# define gsdvt   gsdvt_ 
# define gsdvt2  gsdvt2_
# define gsord   gsord_ 
# define gspos   gspos_ 
# define gsposp  gsposp_ 
# define gsrotm  gsrotm_ 
# define gsvolu  gsvolu_ 
# define gprint  gprint_ 
# define gdinit  gdinit_ 
# define gdopt   gdopt_ 
# define gdraw   gdraw_ 
# define gdrayt  gdrayt_
# define gdrawc  gdrawc_ 
# define gdrawx  gdrawx_ 
# define gdhead  gdhead_ 
# define gdwmn1  gdwmn1_ 
# define gdwmn2  gdwmn2_ 
# define gdwmn3  gdwmn3_ 
# define gdxyz   gdxyz_ 
# define gdman   gdman_ 
# define gdspec  gdspec_ 
# define gdtree  gdtree_ 
# define gdelet  gdelet_ 
# define gdclos  gdclos_ 
# define gdshow  gdshow_ 
# define gdopen  gdopen_ 
# define dzshow  dzshow_ 
# define gsatt   gsatt_ 
# define gfpara  gfpara_
# define gckpar  gckpar_
# define gckmat  gckmat_
# define geditv  geditv_
# define mzdrop  mzdrop_

# define ertrak  ertrak_
# define ertrgo  ertrgo_
 
# define setbomb setbomb_
# define setclip setclip_
# define gcomad gcomad_

# define gbrelm gbrelm_
# define gprelm gprelm_

# define rxgtrak rxgtrak_
# define rxouth  rxouth_
# define rxinh   rxinh_

#else 

# define gzebra  GZEBRA 
# define grfile  GRFILE 
# define gpcxyz  GPCXYZ 
# define ggclos  GGCLOS 
# define glast   GLAST 
# define ginit   GINIT 
# define g3cinit  G3CINIT 
# define gzinit  GZINIT 
# define grun    GRUN 
# define gtrig   GTRIG 
# define gtrigc  GTRIGC 
# define gtrigi  GTRIGI 
# define gfmate  GFMATE 
# define gfpart  GFPART 
# define gftmed  GFTMED 
# define gftmat  GFTMAT
# define gsdk    GSDK 
# define gsmate  GSMATE 
# define gsmixt  GSMIXT 
# define gspart  GSPART 
# define gstmed  GSTMED 
# define gsckov  GSCKOV
# define gstpar  GSTPAR 
# define gfkine  GFKINE 
# define gfvert  GFVERT 
# define gskine  GSKINE 
# define gsvert  GSVERT 
# define gphysi  GPHYSI 
# define gdebug  GDEBUG 
# define gekbin  GEKBIN 
# define gfinds  GFINDS 
# define gsking  GSKING 
# define gskpho  GSKPHO 
# define gsstak  GSSTAK 
# define gsxyz   GSXYZ 
# define gtrack  GTRACK 
# define gtreve  GTREVE 
# define gtreveroot  GTREVEROOT
# define grndm   GRNDM
# define grndmq  GRNDMQ
# define gdtom   GDTOM 
# define glmoth  GLMOTH 
# define gmedia  GMEDIA 
# define gmtod   GMTOD 
# define gsdvn   GSDVN 
# define gsdvn2  GSDVN2 
# define gsdvs   GSDVS 
# define gsdvs2  GSDVS2 
# define gsdvt   GSDVT 
# define gsdvt2  GSDVT2
# define gsord   GSORD 
# define gspos   GSPOS 
# define gsposp  GSPOSP 
# define gsrotm  GSROTM 
# define gsvolu  GSVOLU 
# define gprint  GPRINT 
# define gdinit  GDINIT
# define gdopt   GDOPT 
# define gdraw   GDRAW
# define gdrayt  GDRAYT
# define gdrawc  GDRAWC
# define gdrawx  GDRAWX 
# define gdhead  GDHEAD
# define gdwmn1  GDWMN1
# define gdwmn2  GDWMN2
# define gdwmn3  GDWMN3
# define gdxyz   GDXYZ
# define gdman   GDMAN
# define gdfspc  GDFSPC
# define gdspec  GDSPEC
# define gdtree  GDTREE
# define gdelet  GDELET
# define gdclos  GDCLOS
# define gdshow  GDSHOW
# define gdopen  GDOPEN
# define dzshow  DZSHOW 
# define gsatt   GSATT 
# define gfpara  GFPARA
# define gckmat  GCKMAT
# define geditv  GEDITV
# define mzdrop  MZDROP 

# define ertrak  ERTRAK
# define ertrgo  ERTRGO
 
# define setbomb SETBOMB
# define setclip SETCLIP
# define gcomad  GCOMAD
 
# define gbrelm GBRELM
# define gprelm GPRELM

# define rxgtrak RXGTRAK 
# define rxouth  RXOUTH
# define rxinh   RXINH

#endif 

//
// Geant3 global pointer
//
extern TGeant3 *geant3;

ClassImp(TGeant3f77) 
 
//____________________________________________________________________________ 
TGeant3f77::TGeant3f77() : TGeant3()
{ 
  //
  // Default constructor
  //
} 
 
//____________________________________________________________________________ 
TGeant3f77::TGeant3f77(const char *title, Int_t nwgeant) 
       :TGeant3(title,nwgeant) 
{
} 

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GBASE
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//____________________________________________________________________________ 
extern "C" 
{
  void  g3pvolu(Int_t &); 
  void  g3protm(Int_t &); 
  void  g3ptmed(Int_t &); 
  void  g3pmate(Int_t &); 
  void  g3ppart(Int_t &); 
  void  g3psets(Int_t &); 
  void  g3pvert(Int_t &); 
  void  g3pkine(Int_t &); 
  void  g3pjxyz(Int_t &); 
  void  g3phits(Int_t &); 
  void  g3part(Int_t &); 
  void  g3mate(); 
  void  g3bhsta(); 
  void  g3scank(); 
  void  g3scanu(); 
}

extern "C" {

//____________________________________________________________________________ 
void gfile(const char *filename, const char *option) 
{ 
  //
  //    Routine to open a GEANT/RZ data base. 
  //
  //    LUN logical unit number associated to the file 
  //
  //    CHFILE RZ file name   
  //  
  //    CHOPT is a character string which may be  
  //        N  To create a new file 
  //        U  to open an existing file for update 
  //       " " to open an existing file for read only
  //        Q  The initial allocation (default 1000 records) 
  //           is given in IQUEST(10)
  //        X  Open the file in exchange format
  //        I  Read all data structures from file to memory 
  //        O  Write all data structures from memory to file 
  // 
  // Note:
  //      If options "I"  or "O" all data structures are read or
  //         written from/to file and the file is closed. 
  //      See routine GRMDIR to create subdirectories  
  //      See routines GROUT,GRIN to write,read objects 
  //  
  geant3->Gfile(filename,option);
} 
 
//____________________________________________________________________________ 
void  gpcxyz() 
{ 
  //
  //    Print track and volume parameters at current point
  //
    
    geant3->Gpcxyz(); 
} 
//_____________________________________________________________________________
void  ggclos() 
{ 
  //
  //   Closes off the geometry setting.
  //   Initializes the search list for the contents of each
  //   volume following the order they have been positioned, and
  //   inserting the content '0' when a call to GSNEXT (-1) has
  //   been required by the user.
  //   Performs the development of the JVOLUM structure for all 
  //   volumes with variable parameters, by calling GGDVLP. 
  //   Interprets the user calls to GSORD, through GGORD.
  //   Computes and stores in a bank (next to JVOLUM mother bank)
  //   the number of levels in the geometrical tree and the
  //   maximum number of contents per level, by calling GGNLEV.
  //   Sets status bit for CONCAVE volumes, through GGCAVE.
  //   Completes the JSET structure with the list of volume names 
  //   which identify uniquely a given physical detector, the
  //   list of bit numbers to pack the corresponding volume copy 
  //   numbers, and the generic path(s) in the JVOLUM tree, 
  //   through the routine GHCLOS. 
  //
  geant3->FinishGeometry();
} 
 
//_____________________________________________________________________________
void glast() 
{ 
  //
  // Finish a Geant run
  //
  geant3->Glast(); 
} 
 
//_____________________________________________________________________________
void  gprint(const char *name, Int_t val, const Int_t lname) 
{ 
  //
  // Routine to print data structures
  // CHNAME   name of a data structure
  // 
  if (lname > 0) {
     char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
     if (val == 0) geant3->Gprint(vname); 
     else          geant3->Gprint(vname); 
  } else {
     geant3->Gprint("*");
  }
} 

//_____________________________________________________________________________
void grun() 
{ 
  //
  // Steering function to process one run
  //
  geant3->Grun(); 
} 
 
//_____________________________________________________________________________
void gtrig() 
{ 
  //
  // Steering function to process one event
  //  
  geant3->Gtrig(); 
} 
 
//_____________________________________________________________________________
void gtrigc() 
{ 
  //
  // Clear event partition
  //
  geant3->Gtrigc(); 
} 
 
//_____________________________________________________________________________
void gtrigi() 
{ 
  //
  // Initialises event partition
  //
  geant3->Gtrigi(); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
//_____________________________________________________________________________
void gfmate(Int_t &imat, char *name, Float_t &a, Float_t &z,  
		      Float_t &dens, Float_t &radl, Float_t &absl,
		      Float_t* ubuf, Int_t& nbuf, const Int_t lname) 
{ 
  //
  // Return parameters for material IMAT 
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gfmate(imat, vname,a,z,dens,radl,absl,ubuf,nbuf);
} 
 
//_____________________________________________________________________________
void gfpart(Int_t &ipart, char *name, Int_t &itrtyp,  
		   Float_t &amass, Float_t &charge, Float_t &tlife, const Int_t lname) 
{ 
  //
  // Return parameters for particle of type IPART
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gfpart(ipart, vname, itrtyp, amass, charge, tlife);
} 
 
//_____________________________________________________________________________
void gftmed(Int_t &numed, char *name, Int_t &nmat, Int_t &isvol,  
		   Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd, 
		    Float_t &stemax, Float_t &deemax, Float_t &epsil, 
		    Float_t &stmin, Float_t * /*ubuf*/, Int_t * /*nbuf*/, const Int_t lname) 
{ 
  //
  // Return parameters for tracking medium NUMED
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gftmed(numed, vname, nmat, isvol, ifield, fieldm, tmaxfd, stemax,  
         deemax, epsil, stmin);
}

 
//_____________________________________________________________________________
 void gftmat(Int_t &imate, Int_t &ipart, char *chmeca, Int_t &kdim,
		      Float_t* tkin, Float_t* value, Float_t* pcut,
		      Int_t &ixst)
{ 
  //
  // Return parameters for tracking medium NUMED
  //
  geant3->Gftmat(imate, ipart, chmeca, kdim, tkin, value, pcut, ixst);
} 
 
//_____________________________________________________________________________
void gsdk(Int_t &ipart, Float_t *bratio, Int_t *mode) 
{ 
//  Defines branching ratios and decay modes for standard
//  GEANT particles.
   geant3->Gsdk(ipart,bratio,mode); 
} 
 
//_____________________________________________________________________________
void gsmate(Int_t &imat, const char *name, Float_t &a, Float_t &z,  
		   Float_t &dens, Float_t &radl, Float_t &absl, Float_t */*ubuf*/, Int_t &/*nbuf*/, const Int_t lname) 
{ 
  //
  // Defines a Material
  // 
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gsmate(imat,vname, a, z, dens, radl, absl);
} 
 
//_____________________________________________________________________________
void gsmixt(Int_t &imat, const char *name, Float_t *a, Float_t *z,  
		   Float_t &dens, Int_t &nlmat, Float_t *wmat, const Int_t lname) 
{ 
  //
  //       Defines mixture OR COMPOUND IMAT as composed by 
  //       THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  // 
  //       If NLMAT.GT.0 then WMAT contains the PROPORTION BY
  //       WEIGTHS OF EACH BASIC MATERIAL IN THE MIXTURE. 
  // 
  //       If NLMAT.LT.0 then WMAT contains the number of atoms 
  //       of a given kind into the molecule of the COMPOUND
  //       In this case, WMAT in output is changed to relative
  //       weigths.
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gsmixt(imat,vname, a, z,dens, nlmat,wmat); 
} 
 
//_____________________________________________________________________________
void gspart(Int_t &ipart, const char *name, Int_t &itrtyp,  
		   Float_t &amass, Float_t &charge, Float_t &tlife, Float_t * /*ubuf*/, Int_t & /*nbuf*/, const Int_t lname) 
{ 
  //
  // Store particle parameters
  //
  // ipart           particle code
  // name            particle name
  // itrtyp          transport method (see GEANT manual)
  // amass           mass in GeV/c2
  // charge          charge in electron units
  // tlife           lifetime in seconds
  //
  
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gspart(ipart,vname, itrtyp, amass, charge, tlife);
} 
 
//_____________________________________________________________________________
void gstmed(Int_t &numed, const char *name, Int_t &nmat, Int_t &isvol,  
		      Int_t &ifield, Float_t &fieldm, Float_t &tmaxfd,
		      Float_t &stemax, Float_t &deemax, Float_t &epsil,
		      Float_t &stmin, Float_t * /*ubuf*/, Int_t & /*nbuf*/, const Int_t lname) 
{ 
  //
  //  NTMED  Tracking medium number
  //  NAME   Tracking medium name
  //  NMAT   Material number
  //  ISVOL  Sensitive volume flag
  //  IFIELD Magnetic field
  //  FIELDM Max. field value (Kilogauss)
  //  TMAXFD Max. angle due to field (deg/step)
  //  STEMAX Max. step allowed
  //  DEEMAX Max. fraction of energy lost in a step
  //  EPSIL  Tracking precision (cm)
  //  STMIN  Min. step due to continuos processes (cm)
  //
  //  IFIELD = 0 if no magnetic field; IFIELD = -1 if user decision in GUSWIM;
  //  IFIELD = 1 if tracking performed with GRKUTA; IFIELD = 2 if tracking
  //  performed with G3HELIX; IFIELD = 3 if tracking performed with G3HELX3.
  //  
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gstmed(numed,vname, nmat, isvol, ifield, fieldm, tmaxfd, stemax,
	 deemax, epsil, stmin); 
} 
 
//_____________________________________________________________________________
void gsckov(Int_t &itmed, Int_t &npckov, Float_t *ppckov,
			   Float_t *absco, Float_t *effic, Float_t *rindex)
{ 
  //
  //    Stores the tables for UV photon tracking in medium ITMED 
  //    Please note that it is the user's responsability to 
  //    provide all the coefficients:
  //
  //
  //       ITMED       Tracking medium number
  //       NPCKOV      Number of bins of each table
  //       PPCKOV      Value of photon momentum (in GeV)
  //       ABSCO       Absorbtion coefficients 
  //                   dielectric: absorbtion length in cm
  //                   metals    : absorbtion fraction (0<=x<=1)
  //       EFFIC       Detection efficiency for UV photons 
  //       RINDEX      Refraction index (if=0 metal)
  //
  geant3->Gsckov(itmed,npckov,ppckov,absco,effic,rindex);
}

//_____________________________________________________________________________
void gstpar(Int_t &itmed, const char *name, Float_t &parval, const Int_t lname) 
{ 
  //
  //  To change the value of cut  or mechanism "CHPAR"
  //      to a new value PARVAL  for tracking medium ITMED
  //    The  data   structure  JTMED   contains  the   standard  tracking
  //  parameters (CUTS and flags to control the physics processes)  which
  //  are used  by default  for all  tracking media.   It is  possible to
  //  redefine individually  with GSTPAR  any of  these parameters  for a
  //  given tracking medium. 
  //  ITMED     tracking medium number 
  //  CHPAR     is a character string (variable name) 
  //  PARVAL    must be given as a floating point.
  //
  
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gstpar(itmed,vname, parval); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GCONS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
//_____________________________________________________________________________
void gfkine(Int_t &itra, Float_t *vert, Float_t *pvert, Int_t &ipart,
		      Int_t &nvert, Float_t *ubuf, Int_t &nbuf) 
{ 
  //           Storing/Retrieving Vertex and Track parameters
  //           ---------------------------------------------- 
  //
  //  Stores vertex parameters. 
  //  VERT      array of (x,y,z) position of the vertex 
  //  NTBEAM    beam track number origin of the vertex 
  //            =0 if none exists  
  //  NTTARG    target track number origin of the vertex
  //  UBUF      user array of NUBUF floating point numbers
  //  NUBUF       
  //  NVTX      new vertex number (=0 in case of error). 
  //  Prints vertex parameters.
  //  IVTX      for vertex IVTX.
  //            (For all vertices if IVTX=0) 
  //  Stores long life track parameters.
  //  PLAB      components of momentum 
  //  IPART     type of particle (see GSPART)
  //  NV        vertex number origin of track
  //  UBUF      array of NUBUF floating point user parameters 
  //  NUBUF
  //  NT        track number (if=0 error).
  //  Retrieves long life track parameters.
  //  ITRA      track number for which parameters are requested
  //  VERT      vector origin of the track  
  //  PVERT     4 momentum components at the track origin 
  //  IPART     particle type (=0 if track ITRA does not exist)
  //  NVERT     vertex number origin of the track 
  //  UBUF      user words stored in GSKINE. 
  //  Prints initial track parameters. 
  //  ITRA      for track ITRA 
  //            (For all tracks if ITRA=0) 
  //
  ubuf=0; 
  nbuf=0; 
  geant3->Gfkine(itra,vert,pvert,ipart,nvert); 
} 

//_____________________________________________________________________________
void gfvert(Int_t &nvtx, Float_t *v, Int_t &ntbeam, Int_t &nttarg,
		      Float_t &tofg, Float_t *ubuf, Int_t &nbuf) 
{ 
  //
  //       Retrieves the parameter of a vertex bank
  //       Vertex is generated from tracks NTBEAM NTTARG
  //       NVTX is the new vertex number 
  //
  ubuf=0; 
  nbuf=0; 
  geant3->Gfvert(nvtx,v,ntbeam,nttarg,tofg); 
} 
 
//_____________________________________________________________________________
void gskine(Float_t *plab, Int_t &ipart, Int_t &nv, Float_t *buf,
		      Int_t &nbuf, Int_t &nt) 
{ 
  //
  //       Store kinematics of track NT into data structure
  //       Track is coming from vertex NV
  //
  nt = geant3->Gskine(plab, ipart, nv, buf, nbuf); 
} 
 
//_____________________________________________________________________________
void gsvert(Float_t *v, Int_t &ntbeam, Int_t &nttarg, Float_t *ubuf,
		      Int_t &nbuf, Int_t &nwtx) 
{ 
  //
  //       Creates a new vertex bank 
  //       Vertex is generated from tracks NTBEAM NTTARG 
  //       NVTX is the new vertex number
  //
  nwtx = geant3->Gsvert(v, ntbeam, nttarg, ubuf, nbuf); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GPHYS
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void gphysi() 
{ 
  //
  //       Initialise material constants for all the physics
  //       mechanisms used by GEANT
  //
  geant3->Gphysi(); 
} 
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GTRAK
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 
//_____________________________________________________________________________
void gdebug() 
{ 
  //
  // Debug the current step
  //
  geant3->Gdebug(); 
} 
 
//_____________________________________________________________________________
void gsking(Int_t &igk) 
{ 
  //
  //   Stores in stack JSTAK either the IGKth track of /GCKING/,
  //    or the NGKINE tracks when IGK is 0.
  //
  geant3->Gsking(igk); 
} 
 
//_____________________________________________________________________________
void gskpho(Int_t &igk) 
{ 
  //
  //  Stores in stack JSTAK either the IGKth Cherenkov photon of  
  //  /GCKIN2/, or the NPHOT tracks when IGK is 0.                
  //
  geant3->Gskpho(igk); 
} 
 
//_____________________________________________________________________________
void gsstak(Int_t &iflag) 
{ 
  //
  //   Stores in auxiliary stack JSTAK the particle currently 
  //    described in common /GCKINE/. 
  // 
  //   On request, creates also an entry in structure JKINE :
  //    IFLAG =
  //     0 : No entry in JKINE structure required (user) 
  //     1 : New entry in JVERTX / JKINE structures required (user)
  //    <0 : New entry in JKINE structure at vertex -IFLAG (user)
  //     2 : Entry in JKINE structure exists already (from GTREVE)
  //
  geant3->Gsstak(iflag); 
} 
 
//_____________________________________________________________________________
void gsxyz() 
{ 
  //
  //   Store space point VECT in banks JXYZ 
  //
  geant3->Gsxyz(); 
} 
 
//_____________________________________________________________________________
void gtrack() 
{ 
  //
  //   Controls tracking of current particle 
  //
  geant3->Gtrack(); 
} 
 
//_____________________________________________________________________________
void gtreve() 
{ 
  //
  //   Controls tracking of all particles belonging to the current event
  //
  geant3->Gtreve(); 
} 

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GDRAW
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void gdxyz(Int_t & /*it*/)
{
  //
  // Draw the points stored with Gsxyz relative to track it
  //
  //geant3->Ggdxyz(it);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                        Functions from GGEOM
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void gdtom(Float_t *xd, Float_t *xm, Int_t &iflag) 
{ 
  //
  //  Computes coordinates XM (Master Reference System
  //  knowing the coordinates XD (Detector Ref System)
  //  The local reference system can be initialized by
  //    - the tracking routines and GDTOM used in GUSTEP
  //    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
  //        (inverse routine is GMTOD)
  // 
  //   If IFLAG=1  convert coordinates
  //      IFLAG=2  convert direction cosinus
  //
  geant3->Gdtom(xd, xm, iflag); 
} 
 
//_____________________________________________________________________________
void glmoth(const char* name, Int_t &iunum, Int_t &nlev, Int_t *lvols,
		      Int_t *lindx, const Int_t lname) 
{ 
  //
  //   Loads the top part of the Volume tree in LVOLS (IVO's),
  //   LINDX (IN indices) for a given volume defined through
  //   its name IUDET and number IUNUM.
  // 
  //   The routine stores only upto the last level where JVOLUM
  //   data structure is developed. If there is no development
  //   above the current level, it returns NLEV zero.
  //Int_t *idum=0; 
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Glmoth(vname, iunum, nlev, lvols, lindx); 
} 
 
//_____________________________________________________________________________
void gmtod(Float_t *xm, Float_t *xd, Int_t &iflag) 
{ 
  //
  //       Computes coordinates XD (in DRS) 
  //       from known coordinates XM in MRS 
  //       The local reference system can be initialized by
  //         - the tracking routines and GMTOD used in GUSTEP
  //         - a call to GMEDIA(XM,NUMED)
  //         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
  //             (inverse routine is GDTOM) 
  //
  //        If IFLAG=1  convert coordinates 
  //           IFLAG=2  convert direction cosinus
  //
  geant3->Gmtod(xm, xd, iflag); 
} 
 
//_____________________________________________________________________________
void gsdvn(const char *name, const char *mother, Int_t &ndiv,
		     Int_t &iaxis, const Int_t lname, const Int_t lmother) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  //  NAME   Volume name
  //  MOTHER Mother volume name
  //  NDIV   Number of divisions
  //  IAXIS  Axis value
  //
  //  X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
  //  It divides a previously defined volume.
  //  
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsdvn(vname, vmother, ndiv, iaxis);
} 
 
//_____________________________________________________________________________
void gsdvn2(const char *name, const char *mother, Int_t &ndiv,
		      Int_t &iaxis, Float_t &c0i, Int_t &numed, const Int_t lname, const Int_t lmother) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  // Divides mother into ndiv divisions called name
  // along axis iaxis starting at coordinate value c0.
  // the new volume created will be medium number numed.
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsdvn2(vname, vmother, ndiv, iaxis, c0i, numed);
} 
 
//_____________________________________________________________________________
void gsdvs(const char *name, const char *mother, Float_t &step,
		     Int_t &iaxis, Int_t &numed, const Int_t lname, const Int_t lmother) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsdvs(vname, vmother, step, iaxis, numed);
} 
 
//_____________________________________________________________________________
void gsdvs2(const char *name, const char *mother, Float_t &step,
		      Int_t &iaxis, Float_t &c0, Int_t &numed, const Int_t lname, const Int_t lmother) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsdvs2(vname, vmother, step, iaxis, c0, numed);
} 
 
//_____________________________________________________________________________
void gsdvt(const char *name, const char *mother, Float_t &step,
		     Int_t &iaxis, Int_t &numed, Int_t &ndvmx, const Int_t lname, const Int_t lmother) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  //       Divides MOTHER into divisions called NAME along
  //       axis IAXIS in steps of STEP. If not exactly divisible 
  //       will make as many as possible and will centre them 
  //       with respect to the mother. Divisions will have medium 
  //       number NUMED. If NUMED is 0, NUMED of MOTHER is taken.
  //       NDVMX is the expected maximum number of divisions
  //          (If 0, no protection tests are performed) 
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsdvt(vname, vmother, step, iaxis, numed, ndvmx);
} 

//_____________________________________________________________________________
void gsdvt2(const char *name, const char *mother, Float_t &step,
		      Int_t &iaxis, Float_t &c0, Int_t &numed, Int_t &ndvmx, const Int_t lname, const Int_t lmother) 
{ 
  //
  // Create a new volume by dividing an existing one
  //                                                                    
  //           Divides MOTHER into divisions called NAME along          
  //            axis IAXIS starting at coordinate value C0 with step    
  //            size STEP.                                              
  //           The new volume created will have medium number NUMED.    
  //           If NUMED is 0, NUMED of mother is taken.                 
  //           NDVMX is the expected maximum number of divisions        
  //             (If 0, no protection tests are performed)              
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsdvt2(vname, vmother, step, iaxis, c0, numed, ndvmx);
} 

//_____________________________________________________________________________
void gsord(const char *name, Int_t &iax, const Int_t lname) 
{ 
  //
  //    Flags volume CHNAME whose contents will have to be ordered 
  //    along axis IAX, by setting the search flag to -IAX
  //           IAX = 1    X axis 
  //           IAX = 2    Y axis 
  //           IAX = 3    Z axis 
  //           IAX = 4    Rxy (static ordering only  -> GTMEDI)
  //           IAX = 14   Rxy (also dynamic ordering -> GTNEXT)
  //           IAX = 5    Rxyz (static ordering only -> GTMEDI)
  //           IAX = 15   Rxyz (also dynamic ordering -> GTNEXT)
  //           IAX = 6    PHI   (PHI=0 => X axis)
  //           IAX = 7    THETA (THETA=0 => Z axis)
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gsord(vname, iax);
} 
 
//_____________________________________________________________________________
void gspos(const char *name, Int_t &nr, const char *mother, Float_t &x,
		     Float_t &y, Float_t &z, Int_t &irot, const char *konly, const Int_t lname, const Int_t lmother, const Int_t lkonly) 
{ 
  //
  // Position a volume into an existing one
  //
  //  NAME   Volume name
  //  NUMBER Copy number of the volume
  //  MOTHER Mother volume name
  //  X      X coord. of the volume in mother ref. sys.
  //  Y      Y coord. of the volume in mother ref. sys.
  //  Z      Z coord. of the volume in mother ref. sys.
  //  IROT   Rotation matrix number w.r.t. mother ref. sys.
  //  ONLY   ONLY/MANY flag
  //
  //  It positions a previously defined volume in the mother.
  //  
    
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  char vkonly[24]; strncpy(vkonly,konly,lkonly); vkonly[lkonly] = 0;
  geant3->Gspos(vname, nr, vmother, x, y, z, irot, vkonly);
} 
 
//_____________________________________________________________________________
void gsposp(const char *name, Int_t &nr, const char *mother,  
		      Float_t &x, Float_t &y, Float_t &z, Int_t &irot,
		      const char *konly, Float_t *upar, Int_t &np , const Int_t lname, const Int_t lmother) 
{ 
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vmother[24]; strncpy(vmother,mother,lmother); vmother[lmother] = 0;
  geant3->Gsposp(vname, nr, vmother, x,y, z, irot, konly, upar, np);
}
 
//_____________________________________________________________________________
void gsrotm(Int_t &nmat, Float_t &theta1, Float_t &phi1, Float_t &theta2,
		      Float_t &phi2, Float_t &theta3, Float_t &phi3) 
{ 
  //
  //  nmat   Rotation matrix number
  //  THETA1 Polar angle for axis I
  //  PHI1   Azimuthal angle for axis I
  //  THETA2 Polar angle for axis II
  //  PHI2   Azimuthal angle for axis II
  //  THETA3 Polar angle for axis III
  //  PHI3   Azimuthal angle for axis III
  //
  //  It defines the rotation matrix number IROT.
  //  
  geant3->Gsrotm(nmat, theta1, phi1, theta2, phi2, theta3, phi3); 
} 
 
//_____________________________________________________________________________
void gsvolu(const char *name, const char *shape, Int_t &nmed,  
		      Float_t *upar, Int_t &npar, Int_t &ivolu, const Int_t lname, const Int_t lshape) 
{ 
  //
  //  NAME   Volume name
  //  SHAPE  Volume type
  //  NUMED  Tracking medium number
  //  NPAR   Number of shape parameters
  //  UPAR   Vector containing shape parameters
  //
  //  It creates a new volume in the JVOLUM data structure.
  //  
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vshape[24]; strncpy(vshape,shape,lshape); vshape[lshape] = 0;
  ivolu = geant3->Gsvolu(vname, vshape, nmed, upar, npar);
}
 
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//           T H E    D R A W I N G   P A C K A G E
//           ======================================
//  Drawing functions. These functions allow the visualization in several ways
//  of the volumes defined in the geometrical data structure. It is possible
//  to draw the logical tree of volumes belonging to the detector (DTREE),
//  to show their geometrical specification (DSPEC,DFSPC), to draw them
//  and their cut views (DRAW, DCUT). Moreover, it is possible to execute
//  these commands when the hidden line removal option is activated; in
//  this case, the volumes can be also either translated in the space
//  (SHIFT), or clipped by boolean operation (CVOL). In addition, it is
//  possible to fill the surfaces of the volumes
//  with solid colours when the shading option (SHAD) is activated.
//  Several tools (ZOOM, LENS) have been developed to zoom detailed parts
//  of the detectors or to scan physical events as well.
//  Finally, the command MOVE will allow the rotation, translation and zooming
//  on real time parts of the detectors or tracks and hits of a simulated event.
//  Ray-tracing commands. In case the command (DOPT RAYT ON) is executed,
//  the drawing is performed by the Geant ray-tracing;
//  automatically, the color is assigned according to the tracking medium of each
//  volume and the volumes with a density lower/equal than the air are considered
//  transparent; if the option (USER) is set (ON) (again via the command (DOPT)),
//  the user can set color and visibility for the desired volumes via the command
//  (SATT), as usual, relatively to the attributes (COLO) and (SEEN).
//  The resolution can be set via the command (SATT * FILL VALUE), where (VALUE)
//  is the ratio between the number of pixels drawn and 20 (user coordinates).
//  Parallel view and perspective view are possible (DOPT PROJ PARA/PERS); in the
//  first case, we assume that the first mother volume of the tree is a box with
//  dimensions 10000 X 10000 X 10000 cm and the view point (infinetely far) is
//  5000 cm far from the origin along the Z axis of the user coordinates; in the
//  second case, the distance between the observer and the origin of the world
//  reference system is set in cm by the command (PERSP NAME VALUE); grand-angle
//  or telescopic effects can be achieved changing the scale factors in the command
//  (DRAW). When the final picture does not occupy the full window,
//  mapping the space before tracing can speed up the drawing, but can also
//  produce less precise results; values from 1 to 4 are allowed in the command
//  (DOPT MAPP VALUE), the mapping being more precise for increasing (VALUE); for
//  (VALUE = 0) no mapping is performed (therefore max precision and lowest speed).
//  The command (VALCUT) allows the cutting of the detector by three planes
//  ortogonal to the x,y,z axis. The attribute (LSTY) can be set by the command
//  SATT for any desired volume and can assume values from 0 to 7; it determines
//  the different light processing to be performed for different materials:
//  0 = dark-matt, 1 = bright-matt, 2 = plastic, 3 = ceramic, 4 = rough-metals,
//  5 = shiny-metals, 6 = glass, 7 = mirror. The detector is assumed to be in the
//  dark, the ambient light luminosity is 0.2 for each basic hue (the saturation
//  is 0.9) and the observer is assumed to have a light source (therefore he will
//  produce parallel light in the case of parallel view and point-like-source
//  light in the case of perspective view).
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void gsatt(const char *name, const char *att, Int_t &val, const Int_t lname, const Int_t latt)
{ 
  //
  //  NAME   Volume name
  //  IOPT   Name of the attribute to be set
  //  IVAL   Value to which the attribute is to be set
  //
  //  name= "*" stands for all the volumes.
  //  iopt can be chosen among the following :
  //  
  //     WORK   0=volume name is inactive for the tracking
  //            1=volume name is active for the tracking (default)
  //  
  //     SEEN   0=volume name is invisible
  //            1=volume name is visible (default)
  //           -1=volume invisible with all its descendants in the tree
  //           -2=volume visible but not its descendants in the tree
  //  
  //     LSTY   line style 1,2,3,... (default=1)
  //            LSTY=7 will produce a very precise approximation for
  //            revolution bodies.
  //  
  //     LWID   line width -7,...,1,2,3,..7 (default=1)
  //            LWID<0 will act as abs(LWID) was set for the volume
  //            and for all the levels below it. When SHAD is 'ON', LWID
  //            represent the linewidth of the scan lines filling the surfaces
  //            (whereas the FILL value represent their number). Therefore
  //            tuning this parameter will help to obtain the desired
  //            quality/performance ratio.
  //  
  //     COLO   colour code -166,...,1,2,..166 (default=1)
  //            n=1=black
  //            n=2=red;    n=17+m, m=0,25, increasing luminosity according to 'm';
  //            n=3=green;  n=67+m, m=0,25, increasing luminosity according to 'm';
  //            n=4=blue;   n=117+m, m=0,25, increasing luminosity according to 'm';
  //            n=5=yellow; n=42+m, m=0,25, increasing luminosity according to 'm';
  //            n=6=violet; n=142+m, m=0,25, increasing luminosity according to 'm';
  //            n=7=lightblue; n=92+m, m=0,25, increasing luminosity according to 'm';
  //            colour=n*10+m, m=1,2,...9, will produce the same colour
  //            as 'n', but with increasing luminosity according to 'm';
  //            COLO<0 will act as if abs(COLO) was set for the volume
  //            and for all the levels below it.
  //            When for a volume the attribute FILL is > 1 (and the
  //            option SHAD is on), the ABS of its colour code must be < 8
  //            because an automatic shading of its faces will be
  //            performed.
  //  
  //     FILL  (1992) fill area  -7,...,0,1,...7 (default=0)
  //            when option SHAD is "on" the FILL attribute of any
  //            volume can be set different from 0 (normal drawing);
  //            if it is set to 1, the faces of such volume will be filled
  //            with solid colours; if ABS(FILL) is > 1, then a light
  //            source is placed along the observer line, and the faces of
  //            such volumes will be painted by colours whose luminosity
  //            will depend on the amount of light reflected;
  //            if ABS(FILL) = 1, then it is possible to use all the 166
  //            colours of the colour table, becouse the automatic shading
  //            is not performed;
  //            for increasing values of FILL the drawing will be performed
  //            with higher and higher resolution improving the quality (the
  //            number of scan lines used to fill the faces increases with FILL);
  //            it is possible to set different values of FILL
  //            for different volumes, in order to optimize at the same time
  //            the performance and the quality of the picture;
  //            FILL<0 will act as if abs(FILL) was set for the volume
  //            and for all the levels below it.
  //            This kind of drawing can be saved in 'picture files'
  //            or in view banks.
  //            0=drawing without fill area
  //            1=faces filled with solid colours and resolution = 6
  //            2=lowest resolution (very fast)
  //            3=default resolution
  //            4=.................
  //            5=.................
  //            6=.................
  //            7=max resolution
  //            Finally, if a coloured background is desired, the FILL
  //            attribute for the first volume of the tree must be set
  //            equal to -abs(colo), colo being >0 and <166.
  //  
  //     SET   set number associated to volume name
  //     DET   detector number associated to volume name
  //     DTYP  detector type (1,2)
  //  
//  InitHIGZ();
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vatt[24]; strncpy(vatt,att,latt); vatt[latt] = 0;
  geant3->Gsatt(vname, vatt, val);
} 

//_____________________________________________________________________________
void gfpara(const char *name, Int_t &number, Int_t &intext, Int_t& npar,
			 Int_t& natt, Float_t* par, Float_t* att, const Int_t lname)
{
  //
  // Find the parameters of a volume
  //
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  geant3->Gfpara(vname, number, intext, npar, natt, par, att);
}

//_____________________________________________________________________________
void gckmat(Int_t &itmed, char* natmed)
{
  //
  // Check the parameters of a tracking medium
  //
  geant3->Gckmat(itmed, natmed);
}

//_____________________________________________________________________________
void gdelete(Int_t & /*iview*/)
{ 
  //
  //  IVIEW  View number
  //
  //  It deletes a view bank from memory.
  //
  //geant3->Gdelet(iview);
}
 
//_____________________________________________________________________________
void gdopen(Int_t &iview)
{ 
  //
  //  IVIEW  View number
  //
  //  When a drawing is very complex and requires a long time to be
  //  executed, it can be useful to store it in a view bank: after a
  //  call to DOPEN and the execution of the drawing (nothing will
  //  appear on the screen), and after a necessary call to DCLOSE,
  //  the contents of the bank can be displayed in a very fast way
  //  through a call to DSHOW; therefore, the detector can be easily
  //  zoomed many times in different ways. Please note that the pictures
  //  with solid colours can now be stored in a view bank or in 'PICTURE FILES'
  //
  geant3->Gdopen(iview);
}
 
//_____________________________________________________________________________
void gdclose()
{ 
  //
  //  It closes the currently open view bank; it must be called after the
  //  end of the drawing to be stored.
  //
  //geant3->Gdclos();
}
 
//_____________________________________________________________________________
void gdshow(Int_t & /*iview*/)
{ 
  //
  //  IVIEW  View number
  //
  //  It shows on the screen the contents of a view bank. It
  //  can be called after a view bank has been closed.
  //
  //geant3->Gdshow(iview);
} 

//_____________________________________________________________________________
void gdopt(const char *name,const char *value, const Int_t lname, const Int_t lvalue)
{ 
  //
  //  NAME   Option name
  //  VALUE  Option value
  //
  //  To set/modify the drawing options.
  //     IOPT   IVAL      Action
  //  
  //     THRZ    ON       Draw tracks in R vs Z
  //             OFF (D)  Draw tracks in X,Y,Z
  //             180
  //             360 
  //     PROJ    PARA (D) Parallel projection
  //             PERS     Perspective
  //     TRAK    LINE (D) Trajectory drawn with lines
  //             POIN       " " with markers
  //     HIDE    ON       Hidden line removal using the CG package
  //             OFF (D)  No hidden line removal
  //     SHAD    ON       Fill area and shading of surfaces.
  //             OFF (D)  Normal hidden line removal.
  //     RAYT    ON       Ray-tracing on.
  //             OFF (D)  Ray-tracing off.
  //     EDGE    OFF      Does not draw contours when shad is on.
  //             ON  (D)  Normal shading.
  //     MAPP    1,2,3,4  Mapping before ray-tracing.
  //             0   (D)  No mapping.
  //     USER    ON       User graphics options in the raytracing.
  //             OFF (D)  Automatic graphics options.
  //  
  char vname[24]; strncpy(vname,name,lname); vname[lname] = 0;
  char vvalue[24]; strncpy(vvalue,value,lvalue); vvalue[lvalue] = 0;
  geant3->Gdopt(vname, vvalue);
} 
 
//_____________________________________________________________________________
void gdraw(const char *name,Float_t &theta, Float_t &phi, Float_t &psi,
		    Float_t &u0,Float_t &v0,Float_t &ul,Float_t &vl)
{ 
  //
  //  NAME   Volume name
  //  +
  //  THETA  Viewing angle theta (for 3D projection)
  //  PHI    Viewing angle phi (for 3D projection)
  //  PSI    Viewing angle psi (for 2D rotation)
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  This function will draw the volumes,
  //  selected with their graphical attributes, set by the Gsatt
  //  facility. The drawing may be performed with hidden line removal
  //  and with shading effects according to the value of the options HIDE
  //  and SHAD; if the option SHAD is ON, the contour's edges can be
  //  drawn or not. If the option HIDE is ON, the detector can be
  //  exploded (BOMB), clipped with different shapes (CVOL), and some
  //  of its parts can be shifted from their original
  //  position (SHIFT). When HIDE is ON, if
  //  the drawing requires more than the available memory, the program
  //  will evaluate and display the number of missing words
  //  (so that the user can increase the
  //  size of its ZEBRA store). Finally, at the end of each drawing (with HIDE on),
  //  the program will print messages about the memory used and
  //  statistics on the volumes' visibility.
  //  The following commands will produce the drawing of a green
  //  volume, specified by NAME, without using the hidden line removal
  //  technique, using the hidden line removal technique,
  //  with different linewidth and colour (red), with
  //  solid colour, with shading of surfaces, and without edges.
  //  Finally, some examples are given for the ray-tracing. (A possible
  //  string for the NAME of the volume can be found using the command DTREE).
  //
  geant3->Gdraw(name, theta,phi,psi,u0,v0,ul,vl);
} 
 
//_____________________________________________________________________________
void gdrawc(const char *name,Int_t &axis, Float_t &cut,Float_t &u0,
		     Float_t &v0,Float_t &ul,Float_t &vl)
{ 
  //
  //  NAME   Volume name
  //  CAXIS  Axis value
  //  CUTVAL Cut plane distance from the origin along the axis
  //  +
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  The cut plane is normal to caxis (X,Y,Z), corresponding to iaxis (1,2,3),
  //  and placed at the distance cutval from the origin.
  //  The resulting picture is seen from the the same axis.
  //  When HIDE Mode is ON, it is possible to get the same effect with
  //  the CVOL/BOX function.
  //  
  geant3->Gdrawc(name, axis,cut,u0,v0,ul,vl);
} 
 
//_____________________________________________________________________________
void gdrawx(const char *name,Float_t &cutthe, Float_t &cutphi,
		     Float_t &cutval, Float_t &theta, Float_t &phi, Float_t &u0,
		     Float_t &v0,Float_t &ul,Float_t &vl)
{ 
  //
  //  NAME   Volume name
  //  CUTTHE Theta angle of the line normal to cut plane
  //  CUTPHI Phi angle of the line normal to cut plane
  //  CUTVAL Cut plane distance from the origin along the axis
  //  +
  //  THETA  Viewing angle theta (for 3D projection)
  //  PHI    Viewing angle phi (for 3D projection)
  //  U0     U-coord. (horizontal) of volume origin
  //  V0     V-coord. (vertical) of volume origin
  //  SU     Scale factor for U-coord.
  //  SV     Scale factor for V-coord.
  //
  //  The cut plane is normal to the line given by the cut angles
  //  cutthe and cutphi and placed at the distance cutval from the origin.
  //  The resulting picture is seen from the viewing angles theta,phi.
  //
  geant3->Gdrawx(name, cutthe,cutphi,cutval,theta,phi,u0,v0,ul,vl);
}
 
//_____________________________________________________________________________
void gdhead(Int_t &isel, const char *name, Float_t &chrsiz)
{ 
  //
  //  Parameters
  //  +
  //  ISEL   Option flag  D=111110
  //  NAME   Title
  //  CHRSIZ Character size (cm) of title NAME D=0.6
  //
  //  ISEL =
  //   0      to have only the header lines
  //   xxxxx1 to add the text name centered on top of header
  //   xxxx1x to add global detector name (first volume) on left
  //   xxx1xx to add date on right
  //   xx1xxx to select thick characters for text on top of header
  //   x1xxxx to add the text 'EVENT NR x' on top of header
  //   1xxxxx to add the text 'RUN NR x' on top of header
  //  NOTE that ISEL=x1xxx1 or ISEL=1xxxx1 are illegal choices,
  //  i.e. they generate overwritten text.
  //
  
  geant3->Gdhead(isel,name,chrsiz);
}

//_____________________________________________________________________________
void gdman(Float_t &u, Float_t &v, const char * /*type*/)
{ 
  //
  //  Draw a 2D-man at position (U0,V0)
  //  Parameters
  //  U      U-coord. (horizontal) of the centre of man' R
  //  V      V-coord. (vertical) of the centre of man' R
  //  TYPE   D='MAN' possible values: 'MAN,WM1,WM2,WM3'
  // 
  //   CALL GDMAN(u,v),CALL GDWMN1(u,v),CALL GDWMN2(u,v),CALL GDWMN2(u,v)
  //  It superimposes the picure of a man or of a woman, chosen among
  //  three different ones, with the same scale factors as the detector
  //  in the current drawing.
  //
  
    geant3->Gdman(u,v);
}
 
//_____________________________________________________________________________
void gdspec(const char *name)
{ 
  //
  //  NAME   Volume name
  //
  //  Shows 3 views of the volume (two cut-views and a 3D view), together with
  //  its geometrical specifications. The 3D drawing will
  //  be performed according the current values of the options HIDE and
  //  SHAD and according the current SetClipBox clipping parameters for that
  //  volume.
  //  
  geant3->Gdspec(name);
} 

//_____________________________________________________________________________
void gdtree(const char *name,Int_t &levmax, Int_t &isel)
{ 
  //
  //  NAME   Volume name
  //  LEVMAX Depth level
  //  ISELT  Options
  //
  //  This function draws the logical tree,
  //  Each volume in the tree is represented by a TPaveTree object.
  //  Double-clicking on a TPaveTree draws the specs of the corresponding volume.
  //  Use TPaveTree pop-up menu to select:
  //    - drawing specs
  //    - drawing tree
  //    - drawing tree of parent
  //  
  geant3->Gdtree(name, levmax, isel);
} 

//____________________________________________________________________________
void gekbin()
{
    //
    //=================Create Materials and geometry
    //

   geant3->Gekbin();
}

//____________________________________________________________________________
void ginit()
{
    //
    //=================Create Materials and geometry
    //

   //geant3->Init();
}

//____________________________________________________________________________
void gdinit()
{
    //
    //=================Create Materials and geometry
    //

}

//____________________________________________________________________________
void gpvolu(int &i)
{
   g3pvolu(i);
}
//____________________________________________________________________________
void gprotm(int &i)
{
   g3protm(i);
}
//____________________________________________________________________________
void gpmate(int &i)
{
   g3pmate(i);
}
//____________________________________________________________________________
void gptmed(int &i)
{
   g3ptmed(i);
}
//____________________________________________________________________________
void gppart(int &i)
{
   g3ppart(i);
}
//____________________________________________________________________________
void gpsets(int &i)
{
   g3psets(i);
}
//____________________________________________________________________________
void gpvert(int &i)
{
   g3pvert(i);
}
//____________________________________________________________________________
void gpkine(int &i)
{
   g3pkine(i);
}
//____________________________________________________________________________
void gpjxyz(int &i)
{
   g3pjxyz(i);
}
//____________________________________________________________________________
void gphits(int &i)
{
   g3phits(i);
}
//____________________________________________________________________________
void gpart()
{
   geant3->Gpart();
}
//____________________________________________________________________________
void gmate()
{
   g3mate();
}
//____________________________________________________________________________
void gbhsta()
{
   g3bhsta();
}
//____________________________________________________________________________
void gscank()
{
   g3scank();
}
//____________________________________________________________________________
void gscanu()
{
   g3scanu();
}

//____________________________________________________________________________
void gzinit()
{
    //
    //=================Initialize zebra
    //

   geant3->Gzinit();
}
}
