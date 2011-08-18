//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      a GEANE (sd-system)  track representation
//      (q/p, v',w',v,w,spu) 
//      (v,w) refers to GFDetPlane system
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//      Andrea Fontana       INFN
//
//
//-----------------------------------------------------------

#ifndef GeaneTRACKREP2_HH
#define GeaneTRACKREP2_HH

// Base Class Headers ----------------
#include "GFAbsTrackRep.h"
//#include "FairTrackParP.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op
#include "TVectorT.h"


// Collaborating Class Declarations --
//class FairGeaneProNew;


class GeaneTrackRep2 : public GFAbsTrackRep {
public:

  // Constructors/Destructors ---------
  GeaneTrackRep2();
  GeaneTrackRep2(const GFDetPlane& plane, // will be defined at origin of plane
		 const TVector3& mom,
		 const TVector3& poserr,
		 const TVector3& momerr,
		 int PDGCode);


  virtual ~GeaneTrackRep2();


  virtual GFAbsTrackRep* clone() const {return new GeaneTrackRep2(*this);}
  virtual GFAbsTrackRep* prototype()const{return new GeaneTrackRep2();}

  // Operators
  friend std::ostream& operator<< (std::ostream& s, const GeaneTrackRep2& me);

  // Accessors -----------------------

  // Modifiers

  // Operations ----------------------

  virtual double extrapolate(const GFDetPlane&, TMatrixT<double>& statePred);
  //virtual void extrapolate(const GFDetPlane&, 
  //			   const TMatrixT<double>& stateFrom 
  //			   TMatrixT<double>& stateResult);

  virtual double extrapolate(const GFDetPlane&, 
			   TMatrixT<double>& statePred,
			   TMatrixT<double>& covPred);

  //these two are overriding GFAbsTrackRep methods
  void extrapolateToPoint(const TVector3& pos,
			 TVector3& poca,
			 TVector3& normVec);

  void extrapolateToLine(const TVector3& point1,
	 		 const TVector3& point2,
			 TVector3& poca,
			 TVector3& normVec,
			 TVector3& poca_onwire);


  //  TVector3 getPocaOnLine(const TVector3& p1, 
  //			 const TVector3& p2, 
  //			 bool back=false);

  virtual TVector3 getPos(const GFDetPlane&) ;
  virtual TVector3 getMom(const GFDetPlane&) ;
  virtual void getPosMom(const GFDetPlane&,TVector3& pos,TVector3& mom) ;
  virtual void getPosMomCov(const GFDetPlane& pl,TVector3& pos,TVector3& mom,TMatrixT<double>& cov);
  virtual double getCharge()const {return fCharge;}
  //FairGeaneProNew* getPropagator() {return _geane;}
  int getPDG() {return fPdg;};

  //  void setPropagator(FairGeaneProNew* g){_geane=g;}
  void switchDirection(){;}

  // (-1,0,1) -> (backward prop,decide myself,forward)

private:

  // Private Data Members ------------
  //  FairGeaneProNew* _geane; //!
  //  double _spu; // sign of z-component of momentum

  int fPdg; // pdg code of the particle to be tracked
  int fG3ParticleID; // pdg code of the particle to be tracked
  double fCharge;

  // Private Methods -----------------
  void poca2Line(const TVector3&,const TVector3&,const TVector3&,TVector3&);
  void pocaLine2Line(const TVector3& point1,const TVector3& line1,const TVector3& point2, const TVector3& line2,TVector3& result1,TVector3& result2);

 public:
  ClassDef(GeaneTrackRep2,2)

};


#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
