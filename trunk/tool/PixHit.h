#ifndef PIXHIT_HH
#define PIXHIT_HH

#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"

typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class PixHit : public PlanarRecoHit {
public:

  PixHit();
  PixHit(const TVector3& point,const TVector3& norm, const TVector3& u, double res, bool smear=false);

  virtual ~PixHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

private:
  static const int NparHitRep = 2;

public:
  ClassDef(PixHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
