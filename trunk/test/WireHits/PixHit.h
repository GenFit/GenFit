
#ifndef PIXHIT_HH
#define PIXHIT_HH

// Base Class Headers ----------------
#include "GFRecoHitIfc.h"
#include "GFPlanarHitPolicy.h"
#include "GFDetPlane.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op


// Collaborating Class Declarations --

typedef GFRecoHitIfc<GFPlanarHitPolicy> PlanarRecoHit;

class PixHit : public PlanarRecoHit {
public:

  // Constructors/Destructors ---------
  PixHit();
  PixHit(TVector3 point,double res);
  PixHit(const GFDetPlane& pl,double res);

  virtual ~PixHit();

  virtual GFAbsRecoHit* clone();
  
  // Operations ----------------------
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector);


private:
  static const int NparHitRep = 2;

  // Private Methods -----------------

public:
  ClassDef(PixHit,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
