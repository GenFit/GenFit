#ifndef POINTHIT_HH
#define POINTHIT_HH

// Base Class Headers ----------------
#include "GFRecoHitIfc.h"
#include "GFSpacepointHitPolicy.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --

typedef GFRecoHitIfc<GFSpacepointHitPolicy> SpacepointRecoHit;

class PointHit : public SpacepointRecoHit {
public:

  // Constructors/Destructors ---------
  PointHit();
  PointHit(const TVector3&);

  virtual ~PointHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);


private:

  // Private Data Members ------------
  //static const int NparHitRep = 3;

  // Private Methods -----------------

public:
  ClassDef(PointHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
