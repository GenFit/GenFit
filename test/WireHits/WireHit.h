#ifndef WIREHIT_HH
#define WIREHIT_HH

// Base Class Headers ----------------
#include "GFRecoHitIfc.h"
#include "GFWireHitPolicy.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --

typedef GFRecoHitIfc<GFWireHitPolicy> WireRecoHit;

class WireHit : public WireRecoHit {
public:

  // Constructors/Destructors ---------
  WireHit();
  WireHit(const TVector3& p1,const TVector3& p2,double driftLength,double res);

  virtual ~WireHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);


private:

  // Private Data Members ------------
  static const int NparHitRep = 7;

  // Private Methods -----------------

public:
  ClassDef(WireHit,1)

};

#endif

