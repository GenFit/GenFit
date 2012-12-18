#ifndef WIREPOINTHIT_HH
#define WIREPOINTHIT_HH

#include "RecoHits/GFAbsWirepointHit.h"


class WirePointHit : public GFAbsWirepointHit {
 public:

  WirePointHit();
  WirePointHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double z,
               double resDrift, double resZ, bool smear=false);

  virtual ~WirePointHit();

  virtual GFAbsRecoHit* clone();
  
  virtual const TMatrixT<double>& getHMatrix(const GFAbsTrackRep* rep);

 public:
  ClassDef(WirePointHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
