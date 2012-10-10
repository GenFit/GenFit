#ifndef PseudoSpacePointWireHit_HH
#define PseudoSpacePointWireHit_HH

#include "GFRecoHitIfc.h"
#include "GFPseudoSpacepointWireHitPolicy.h"

typedef GFRecoHitIfc<GFPseudoSpacepointWireHitPolicy> PseudoSpacePointWireRecoHit;

class PseudoSpacePointWireHit : public PseudoSpacePointWireRecoHit {
public:

  PseudoSpacePointWireHit();
  PseudoSpacePointWireHit(const TVector3& pos, const TVector3& wireDir,
                          double resPerp, double resWire, bool smear=false);

  virtual ~PseudoSpacePointWireHit();

  virtual GFAbsRecoHit* clone();
  
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

  void setWireDirection(TVector3 &v){
    fPolicy.setWireDirection(v);
  }

private:
  static const int NparHitRep = 3;

public:
  ClassDef(PseudoSpacePointWireHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
