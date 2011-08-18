
// Configuration macro for Geane VirtualMC 
#include<iostream>
void Geane()
{

  GeaneMCApplication* app = new GeaneMCApplication();

  //the following gMC3 config code was taken from the Pandaroot framework
  //and was originally written bu M. Al-Turany (GSI, Darmstadt)
  gSystem->Load("libgeant321");
  gMC3= new  TGeant3TGeo("C++ Interface to Geant3");
  // ******* GEANEconfiguration for simulated Runs  *******
  gMC3->SetDEBU(0, 0, 1);
  gMC3->SetSWIT(4, 10);  
  gMC3->SetDCAY(0);
  gMC3->SetPAIR(0);
  gMC3->SetCOMP(0);
  gMC3->SetPHOT(0);
  gMC3->SetPFIS(0); 
  gMC3->SetDRAY(0);
  gMC3->SetANNI(0);
  gMC3->SetBREM(1);
  gMC3->SetMUNU(0);
  gMC3->SetCKOV(0);
  gMC3->SetHADR(0);
  gMC3->SetLOSS(4);
  gMC3->SetMULS(1); 	    //1=Moliere,3=Gaussian
  gMC3->SetRAYL(0);
  gMC3->SetSTRA(0);
  gMC3->SetAUTO(1);         //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
  gMC3->SetABAN(0);         //Restore 3.16 behaviour for abandoned tracks
  gMC3->SetOPTI(0);         //Select optimisation level for GEANT geometry searches (0,1,2)
  gMC3->SetERAN(5.e-7);
  Float_t cut =  1.e-3;                               // 1 MeV cut by default
  Float_t cutd = 1.e4 ;                               // 10 TeV - Threshold for delta-rays
  Float_t cutb = cutd;                                // 10 TeV - Cut for bremsstrahlung
  Float_t tofmax = 1.e10;                             // seconds
  Float_t usrcuts[5] = {0.,0.,0.,0.,0.};              // usercuts
  Float_t gcalpha = 0.999;                            // Optimal value for alpha 
  // set cuts here 
  //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
  gMC3->SetCUTS(cut,  		// CUTGAM = gammas
	 	cut,   	        // CUTELE = electrons
		cut,   	        // CUTNEU = neutral hadrons
		cut,   	        // CUTHAD = charged hadrons
		cut,   	        // CUTMUO = muons
		cutb,  		// BCUTE  = electron bremsstrahlung
		cutb,  		// BCUTM  = muon bremsstrahlung
		cutd,  		// DCUTE  = delta rays by electrons
		cutd,  		// DCUTM  = delta rays by muons
		cutb,   	// PPCUTM = pair production by muons
		tofmax, 	// TOFMAX = time of flight cut
		usrcuts);   
  gMC3->SetECut(gcalpha);
  //end of config from Pandroot

  app->InitMC();

}

