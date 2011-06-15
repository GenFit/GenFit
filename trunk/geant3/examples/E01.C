// $Id: E01.C 220 2007-11-19 16:08:06Z rdm $
//
// Macro for running Example01 with Geant3 
// Before running this macro, the libexampl01.so library
// must have been built. To build it, go to your geant4_vmc/examples directory
// and run make.
// Note that this macro is a simplified version of the equivalent macro
// in the geant4_vmc/examples/E01 directory

{
  // Load basic libraries
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libEG"); 
  gSystem->Load("libEGPythia6");
  gSystem->Load("libPythia6");  
  new TGeoManager("E01","test"); 
  
  // Load Geant3 libraries
  gSystem->Load("../lib/tgt_linux/libdummies");
  gSystem->Load("../lib/tgt_linux/libgeant321");
  
  // Load this example library
  gSystem->Load("~/geant4_vmc/lib/tgt_linux/libexample01");

  // MC application
  Ex01MCApplication* appl 
    = new Ex01MCApplication("Example01", "The example01 VMC application");

  appl->InitMC("E01_g3Config.C");
  appl->RunMC(1);
}  
