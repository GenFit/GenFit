#ifndef WIREPOINTHIT_HH
#define WIREPOINTHIT_HH

#include "GFRecoHitIfc.h"
#include "GFWirepointHitPolicy.h"

typedef GFRecoHitIfc<GFWirepointHitPolicy> WirePointRecoHit;

class WirePointHit : public WirePointRecoHit {
public:

  WirePointHit();
  WirePointHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double z,
               double resDrift, double resZ, bool smear=false);

  virtual ~WirePointHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

private:
  static const int NparHitRep = 8;

public:
  ClassDef(WirePointHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
