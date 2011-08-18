#ifndef SLTRACKREP_HH
#define SLTRACKREP_HH
#include "GFAbsTrackRep.h"
class SlTrackRep : public GFAbsTrackRep {
public:

  // Constructors/Destructors ---------
  SlTrackRep();
  SlTrackRep(const TMatrixT<double>&, const TMatrixT<double>&);
  SlTrackRep(const TMatrixT<double>&, const TMatrixT<double>&,const double);
  SlTrackRep(const GFDetPlane&,const TMatrixT<double>&, const TMatrixT<double>&);
  SlTrackRep(const TVector3& pos, const TVector3& dir);


  virtual ~SlTrackRep();


  virtual GFAbsTrackRep* clone() const {return new SlTrackRep(*this);}
  virtual GFAbsTrackRep* prototype()const{return new SlTrackRep();}

  void setReferencePlane(const GFDetPlane& pl) {fRefPlane=pl;}

  // Operations ----------------------

  virtual double extrapolate(const GFDetPlane&, 
			   TMatrixT<double>& statePred,
			   TMatrixT<double>& covPred);
  virtual double extrapolate(const GFDetPlane&, 
			   TMatrixT<double>& statePred);


  void extrapolateToPoint(const TVector3& pos,
			 TVector3& poca,
			 TVector3& dirInPoca);

  void extrapolateToLine(const TVector3& point1,
	 		 const TVector3& point2,
			 TVector3& poca,
			 TVector3& dirInPoca,
			 TVector3& poca_onwire);


  virtual TVector3 getPos(const GFDetPlane&) ;
  virtual TVector3 getMom(const GFDetPlane&);

  virtual void getPosMom(const GFDetPlane&,TVector3& pos,TVector3& mom) ;
  virtual double getCharge()const {return 0;}

  void switchDirection(){_backw=-_backw;}

private:

  // Private Data Members ------------
  int _backw; // (-1,0,1) -> (backward prop,decide myself,forward)

 public:
  ClassDef(SlTrackRep,1);

};


#endif

