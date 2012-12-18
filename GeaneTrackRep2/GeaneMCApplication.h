
#ifndef GEANE_MC_APPLICATION_H
#define GEANE_MC_APPLICATION_H

#include "TVirtualMCApplication.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include <iostream>

//#include "AbsBField.h"

class GeaneMCApplication : public TVirtualMCApplication
{

  public:
    GeaneMCApplication();
    virtual ~GeaneMCApplication(){;}
    void InitMC();
    /** Construct user geometry */


    virtual void          Field(const Double_t* x, Double_t* b) const;

    //Geant action methods need to be overridden
    virtual void          ConstructGeometry();
    virtual void          FinishEvent(){;}
    virtual void          FinishPrimary(){;}
    virtual void          FinishRun(){;}
    virtual void          GeneratePrimaries(){;}
    virtual void          InitGeometry(){;}
    virtual void          PostTrack(){;}
    virtual void          PreTrack(){;}
    virtual void          BeginEvent(){;}
    virtual void          BeginPrimary(){;}
    virtual void          Stepping(){;}
private:
public:
    ClassDef(GeaneMCApplication,1)
};



#endif
