// $Id: TG3Application.h 220 2007-11-19 16:08:06Z rdm $
//
// Class TG3Application
// ----------------------- 
// Implementation of the TVirtualMCApplication
//
// by Rene Brun 30/12/2002

#ifndef TG3Application_h
#define TG3Application_h

#include <TVirtualMCApplication.h>
#include <TVirtualMCStack.h>

class TG3Application : public TVirtualMCApplication
{
  public:
    TG3Application(const char *name, const char *title);
    TG3Application();
    virtual ~TG3Application();
  
    // static access method
    static TG3Application* Instance(); 

    // methods
    void InitMC(const char *setup);
    void RunMC(Int_t nofEvents);
    void FinishRun();
 
    virtual void ConstructGeometry();
    virtual void InitGeometry();
    virtual void GeneratePrimaries();
    virtual void BeginEvent();
    virtual void BeginPrimary();
    virtual void PreTrack();
    virtual void Stepping();
    virtual void PostTrack();
    virtual void FinishPrimary();
    virtual void FinishEvent();
    
    virtual Double_t TrackingRmax() const;
    virtual Double_t TrackingZmax() const;
    virtual void Field(const Double_t* x, Double_t* b) const;

  private:
    // methods
    void ConstructMaterials();
    void ConstructVolumes();
  
    // data members
    TVirtualMCStack*  fStack;

  ClassDef(TG3Application,1)  //dummy Interface to G3 MonteCarlo application
};

// inline functions

inline TG3Application* TG3Application::Instance()
{ return (TG3Application*)(TVirtualMCApplication::Instance()); }

#endif //TG3Application_h

