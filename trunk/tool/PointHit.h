#ifndef POINTHIT_HH
#define POINTHIT_HH

#include "GFRecoHitIfc.h"
#include "GFSpacepointHitPolicy.h"

typedef GFRecoHitIfc<GFSpacepointHitPolicy> SpacepointRecoHit;

class PointHit : public SpacepointRecoHit {
public:

  // Constructors/Destructors ---------
  PointHit();
  PointHit(const TVector3& point,const double& res);
  PointHit(const TVector3& point,const TVector3& res);

  virtual ~PointHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);

private:
  static const int NparHitRep = 3;
public:
  ClassDef(PointHit,1)
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
