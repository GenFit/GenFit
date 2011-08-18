#ifndef STRIPHIT_HH
#define STRIPHIT_HH

#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"

typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class StripHit : public PlanarRecoHit {
public:

  StripHit();
  StripHit(const TVector3& point,const TVector3& norm, const TVector3& u,double res,double smear,int smearFlag);

  virtual ~StripHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);



public:
  ClassDef(StripHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
