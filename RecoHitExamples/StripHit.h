#ifndef STRIPHIT_HH
#define STRIPHIT_HH

#include "RecoHits/GFAbsPlanarHit.h"


class StripHit : public GFAbsPlanarHit {
public:

  StripHit();
  StripHit(const TVector3& point,const TVector3& norm, const TVector3& u, double res, bool smear=false);

  virtual ~StripHit();

  virtual GFAbsRecoHit* clone();
  
  virtual const TMatrixT<double>& getHMatrix(const GFAbsTrackRep* rep);

private:
  static const int NparHitRep = 1;

public:
  ClassDef(StripHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
