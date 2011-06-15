// $Id: TG3Application.cxx 220 2007-11-19 16:08:06Z rdm $
//
// Class TG3Application
// ----------------------- 
// Implementation of the TVirtualMCApplication
//
// by Rene Brun 30/12/2002

#include "TG3Application.h"
#include "TGeant3f77.h"
//#include "TG3Stack.h"

#include <TROOT.h>
#include <TInterpreter.h>
#include <TVirtualMC.h>
#include <TLorentzVector.h>

#include <Riostream.h>

ClassImp(TG3Application)

//_____________________________________________________________________________
TG3Application::TG3Application(const char *name, const char *title) 
  : TVirtualMCApplication(name,title),
    fStack(0)
{
//
// Standard constructor
//

  new TGeant3f77("g",0);
  
  // create a user stack
  //fStack = new TG3Stack(100);  
}

//_____________________________________________________________________________
TG3Application::TG3Application()
  : TVirtualMCApplication(),
    fStack(0)
{    
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
TG3Application::~TG3Application() 
{
  //
  // Destructor  
  //
  
  //delete fStack;
}

//
// private
//

//_____________________________________________________________________________
void TG3Application::ConstructMaterials()
{
  //
  // Materials
  //

}


//_____________________________________________________________________________
void TG3Application::ConstructVolumes()
{

}

//
// public
//

//_____________________________________________________________________________
void TG3Application::InitMC(const char* /*setup*/)
{    
  //
  // Initialize MC.
  //

  gMC->Init();
  gMC->BuildPhysics();  
}

//_____________________________________________________________________________
void TG3Application::RunMC(Int_t /*nofEvents*/)
{    
  //
  // MC run.
  //


  //gMC->ProcessRun(nofEvents);
  //FinishRun();
}

//_____________________________________________________________________________
void TG3Application::FinishRun()
{    
  //
  // Finish MC run.
  //

 // UGLAST
}

//_____________________________________________________________________________
void TG3Application::ConstructGeometry()
{    
  //
  // Construct geometry using TVirtualMC functions.
  //

  //ConstructMaterials();  
  //ConstructVolumes();  
}

//_____________________________________________________________________________
void TG3Application::InitGeometry()
{    
  //
  // Initialize geometry
  //
  
  // Nothing needed in this example
}

//_____________________________________________________________________________
void TG3Application::GeneratePrimaries()
{    
  //
  // Fill the user stack (derived from TVirtualMCStack) with primary particles.
  //
  
}

//_____________________________________________________________________________
void TG3Application::BeginEvent()
{    
  //
  // User actions at beginning of event
  //

  // nothing to be done this example
}

//_____________________________________________________________________________
void TG3Application::BeginPrimary()
{    
  //
  // User actions at beginning of a primary track
  //

  // nothing to be done this example
}

//_____________________________________________________________________________
void TG3Application::PreTrack()
{    
  //
  // User actions at beginning of each track
  //
}

//_____________________________________________________________________________
void TG3Application::Stepping()
{    
  //
  // User actions at each step
  //
  
}

//_____________________________________________________________________________
void TG3Application::PostTrack()
{    
  //
  // User actions after finishing of each track
  //

  // nothing to be done this example
}

//_____________________________________________________________________________
void TG3Application::FinishPrimary()
{    
  //
  // User actions after finishing of a primary track
  //

  // nothing to be done this example
}

//_____________________________________________________________________________
void TG3Application::FinishEvent()
{    
  //
  // User actions after finishing of an event
  //

  // nothing to be done this example
} 

//_____________________________________________________________________________
Double_t TG3Application::TrackingRmax() const
{ 
  //
  // No limit
  //

  return DBL_MAX; 
}

//_____________________________________________________________________________
Double_t TG3Application::TrackingZmax() const
{ 
  //
  // No limit
  //

  return DBL_MAX; 
}

//_____________________________________________________________________________
void TG3Application::Field(const Double_t* /* x */, Double_t* b) const
{
  // 
  // No magnetic field.
  //
  
   b[0] = 0.;
   b[1] = 0.;
   b[2] = 0.;
}
