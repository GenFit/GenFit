#ifndef SpacepointHit_HH
#define SpacepointHit_HH

#include "RecoHits/GFAbsSpacepointHit.h"

class SpacepointHit : public GFAbsSpacepointHit {
public:

  // Constructors/Destructors ---------
  SpacepointHit();
  SpacepointHit(const TVector3& point, const double& res, bool smear=false);
  SpacepointHit(const TVector3& point, const TVector3& res, bool smear=false);

  virtual ~SpacepointHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual const TMatrixT<double>& getHMatrix(const GFAbsTrackRep* rep);

public:
  ClassDef(SpacepointHit,1)
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
