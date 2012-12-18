#ifndef WIREHIT_HH
#define WIREHIT_HH

#include "RecoHits/GFAbsWireHit.h"


class WireHit : public GFAbsWireHit {
public:

  WireHit();
  WireHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double res, bool smear=false);

  virtual ~WireHit();

  virtual GFAbsRecoHit* clone();
  
  virtual const TMatrixT<double>& getHMatrix(const GFAbsTrackRep* rep);

public:
  ClassDef(WireHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
