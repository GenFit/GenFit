#ifndef WIREHIT_HH
#define WIREHIT_HH

#include "GFRecoHitIfc.h"
#include "GFWireHitPolicy.h"

typedef GFRecoHitIfc<GFWireHitPolicy> WireRecoHit;

class WireHit : public WireRecoHit {
public:

  WireHit();
  WireHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double res, bool smear=false);

  virtual ~WireHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

private:
  static const int NparHitRep = 7;

public:
  ClassDef(WireHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
