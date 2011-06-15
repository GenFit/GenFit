// $Id: E03.C 220 2007-11-19 16:08:06Z rdm $
//
// Macro for running Example03 with Geant3 
// Before running this macro, the libexampl03.so library
// must have been built. To build it, go to your geant4_vmc/examples directory
// and run make.
// Note that this macro is a simplified version of the equivalent macro
// in the geant4_vmc/examples/E03 directory

{
  // Load basic libraries
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libEG"); 
  gSystem->Load("libEGPythia6");
  gSystem->Load("libPythia6");  
  new TGeoManager("E02","test"); 
  
  // Load Geant3 libraries
  gSystem->Load("../lib/tgt_linux/libdummies");
  gSystem->Load("../lib/tgt_linux/libgeant321.so");
  
  // Load this example library
  gSystem->Load("~/geant4_vmc/lib/tgt_linux/libexample03");

  // MC application
  Ex03MCApplication* appl 
    =  new Ex03MCApplication("Example03", "The example03 VMC application");
  appl->GetPrimaryGenerator()->SetNofPrimaries(20);
  appl->SetPrintModulo(1);

  appl->InitMC("E03_g3Config.C");

  // visualization setting
  gMC->Gsatt("*", "seen", 16);
  gMC->Gsatt("ABSO", "seen", 5);
  gMC->Gsatt("GAPX", "seen", 6);

  appl->RunMC(50);
}  
