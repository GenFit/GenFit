/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

/** @addtogroup genfit
 * @{
 */

#ifndef GFABSTRACKREP_H
#define GFABSTRACKREP_H

#include <vector>
#include <list>
#include <iostream>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <Math/ProbFunc.h>

#include "GFDetPlane.h"

class TVector3;
class GFAbsRecoHit;    

/** @brief Base Class for genfit track representations. 
 * Defines interface for track parameterizations.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * It is important to understand the difference between a track and a 
 * track representation in genfit:
 * - A track representation is a specific parameterization of a trajectory. 
 * It contains the parameters that describe the track at some point and 
 * code for the extrapolation of the track parameters through space. 
 * - A Track is a collection of RecoHits (see GFAbsRecoHit) plus a collection
 * of track representation objects. The hits can be from different detectors.
 * There can be several representations of the same track. This makes it 
 * possible to perform several fits in parallel, for example to compare 
 * different parameterizations or to fit different particle hypotheses.
 *
 * All track representations must inherit GFAbsTrackRep to be available
 * in genfit. Algorithms in genfit use this class as interface to 
 * access track parameters
 *
 * Provides:
 *  - Matrix objects to store track parameters
 *  - ... and covariances
 *  - interface to track extrapolation code
 *
 * The track extrapolation engine can be exchanged in genfit. 
 * Or one can even use more than one engine in parallel!
 * In order to use a track extrapolation engine (like e.g. GEANE) with genfit
 * one has to write a TrackRep class that inherits from GFAbsTrackRep. This makes
 * it possible to uses different track extrapolation codes within a unified
 * framework without major changes in the detector code.
 *
 * There is only one thing one has to do to use a specific track 
 * representation together with the hits from a detector: 
 * add the respective code in the GFAbsRecoHit::getHMatrix method implementation
 * of the RecoHit in question.
 */

class GFAbsTrackRep : public TObject{

 public:
  virtual GFAbsTrackRep* clone() const = 0;

  virtual GFAbsTrackRep* prototype() const = 0;
  
  //! returns the track-length spanned in this extrapolation
  /*! There is a default implementation in GFAbsTrackRep.cxx which just drops
     the predicted covariance. If your trackrep has a way to extrapolate
     without giving a correct cov (that would be faster probably), please
     overwrite it.
     This method does NOT alter the state of the object!
  */
  virtual double extrapolate(const GFDetPlane& plane, TVectorD& statePred);


 public:

  GFAbsTrackRep();
  GFAbsTrackRep(int);
  virtual ~GFAbsTrackRep();


  //! This method is to extrapolate the track to point of closest approach to a point in space
  /*! There is an empty implementation of this method in GFAbsTrackRep.cxx,
      which will just abort with an error message. One can overwrite this
      method if one wishes to implement a track representation, which should
      have this feature. An example of an experiment in which you would not
      need this feature would be track fitting (not so much vertexing) in
      an experiment with only planar trackers like silicons or planar wire
      chambers and such. An example where you would need it, would be a TPC
      where you have to fit the track to space points, or other drift chambers
      with complicated hit topology.
      This method does NOT alter the state of the object!
   */
  virtual double extrapolateToPoint(const TVector3& point,
         TVector3& poca,
         TVector3& normVec);

  //! 
  /** @brief This method extrapolates to the point of closest approach to a line
   * 
   * This method extrapolates to the POCA to a line, i.e. a wire. There
   * is a default implementation just like for the extrapolateToPoca for
   * trackReps which do not need this feature, which will abort the
   * execution if it is ever called.
   * This method does NOT alter the state of the object!
   */
  virtual double extrapolateToLine(const TVector3& point1,
         const TVector3& point2,
         TVector3& poca,
         TVector3& normVec,
         TVector3& poca_onwire);
  
  
  //! make step of h cm along the track
  /*! There is an empty implementation in GFAbsTrackRep.cxx which will abort
      (see one of the extrapolate methods above). This can be overwritten,
      if this feature is needed.
      This method does NOT alter the state of the object!
  */
  virtual double stepalong(double h,
           TVector3& point,
           TVector3& dir);

  //! Extrapolates the track to the given detector-plane
  /*! Results are put into statePred and covPred
      This method does NOT alter the state of the object!
   */ 
  virtual double extrapolate(const GFDetPlane& plane,
           TVectorD& statePred,
           TMatrixDSym& covPred)=0;
  
  //! This changes the state and cov and plane of the rep
  /*! This method extrapolates to to the plane and sets the results of state,
      cov and also plane in itself.
  */
  double extrapolate(const GFDetPlane& plane);

  //! returns dimension of state vector
  unsigned int getDim() const {return fDimension;}  
  
  virtual void Print(const Option_t* = "") const;

  inline const TVectorD& getState() const {
    return fState;
  }
  inline const TMatrixDSym& getCov() const {
    return fCov;
  }

  double getStateElem(int i) const {return fState(i);}
  double getCovElem(int i, int j) const {return fCov(i,j);}

  

  virtual TVector3 getPos(const GFDetPlane& pl)=0; 
  virtual TVector3 getMom(const GFDetPlane& pl)=0; 
  virtual void getPosMom(const GFDetPlane& pl,TVector3& pos,TVector3& mom)=0;

  //! method which gets position, momentum and 6x6 covariance matrix
  /*! 
   * default implementation in cxx file, if a ConcreteTrackRep can 
   * not implement this functionality
   */
  virtual void getPosMomCov(const GFDetPlane& pl, TVector3& pos, TVector3& mom, TMatrixDSym& cov);

  virtual double getCharge()const =0;

  //! method which gets pdg particle code
  /*!
   * default implementation in cxx; should be implemented in a ConcreteTrackRep
   */
  virtual int getPDG();

  TVector3 getPos() {return getPos(fRefPlane);}
  TVector3 getMom() {return getMom(fRefPlane);}

  void getPosMomCov(TVector3& pos, TVector3& mom, TMatrixDSym& c){
    getPosMomCov(fRefPlane,pos,mom,c);
  }

  inline const TVectorD& getFirstState() const {
    return fFirstState;
  }
  inline const TMatrixDSym& getFirstCov() const {
    return fFirstCov;
  }
  inline GFDetPlane getFirstPlane() const {
    return fFirstPlane;
  }
  inline const TVectorD& getLastState() const {
    return fLastState;
  }
  inline const TMatrixDSym& getLastCov() const {
    return fLastCov;
  }
  inline GFDetPlane getLastPlane() const {
    return fLastPlane;
  }
  inline double getChiSqu() const {
    return fChiSqu;
  }
  inline double getForwardChiSqu() const {
    return fForwardChiSqu;
  }
  //! returns chi2/ndf
  inline double getRedChiSqu() const {
    if(getNDF()>0) {
      return getChiSqu()/getNDF();
    } else {
      return 0;
    }
  }
  /** @brief returns the p value
   */
  inline double getPVal() const {
    return ROOT::Math::chisquared_cdf_c(getChiSqu(), getNDF()); //Prob only works with integer ndfs so here chisquared_cdf_c is needed
  }
  inline double getNDF() const {
    if(fNdf>getDim()){
      return fNdf-getDim();
    } else {
      return 0;
    }
  }
  /** @brief  X/X0 (total fraction of radiation length passed), cumulated during last extrapolation.
   *  The fitter has to reset XX0 via resetXX0()
   */
  inline double getXX0() const {
    return fXX0;
  }
  /** @brief Puts the track representation in a given state.
   * 
   * This is used to update the track representation after the update of the 
   * Kalman filter was calculated.\n
   * \n
   * IMPORTANT: One should be able to set the track representation to arbitrary
   * values using this method. If the track representation needs additional 
   * information beside the state vector, the plane and the covariance, it has
   * to be handed over via the "aux" Matrix. GFAbsTrackRep::getAuxInfo() should
   * return the appropriate information. This is mandatory if smoothing is used.
   */
  virtual void setData(const TVectorD& st, const GFDetPlane& pl, const TMatrixDSym* cov=NULL, const TMatrixD* aux=NULL){
    fState = st;
    if (&fRefPlane != &pl) fRefPlane = pl; // Copy reference plane only if it changed.
    if(cov!=NULL) fCov = *cov;
    static_cast<void>(aux);
  }
  inline void setCov(const TMatrixDSym& aCov) {
    fCov=aCov;
  }
  inline void setFirstState(const TVectorD& aState) {
    fFirstState = aState;
  }
  inline void setFirstCov(const TMatrixDSym& aCov) {
    fFirstCov = aCov;
  }
  inline void setFirstPlane(const GFDetPlane& aPlane) {
    fFirstPlane = aPlane;
  }
  inline void setLastState(const TVectorD& aState) {
    fLastState = aState;
  }
  inline void setLastCov(const TMatrixDSym& aCov) {
    fLastCov = aCov;
  }
  inline void setLastPlane(const GFDetPlane& aPlane) {
    fLastPlane = aPlane;;
  }

  //! method which sets position, momentum and 6x6 covariance matrix
  /*!
   * default implementation in cxx file, if a ConcreteTrackRep can
   * not implement this functionality
   */
  virtual void setPosMomCov(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov);

  const GFDetPlane& getReferencePlane() const {return fRefPlane;}

  inline void setChiSqu(double aChiSqu) {
    fChiSqu = aChiSqu;
  }
  inline void setForwardChiSqu(double aChiSqu) {
    fForwardChiSqu = aChiSqu;
  }
  inline void setNDF(double n) {
    fNdf = n;
  }
  inline void addChiSqu(double aChiSqu) {
    fChiSqu += aChiSqu;
  }
  inline void addForwardChiSqu(double aChiSqu) {
    fForwardChiSqu += aChiSqu;
  }
  inline void addNDF(double n) {
    fNdf += n;
  }
  inline void setStatusFlag(int _val) {
    fStatusFlag = _val;
  }
  
  virtual void switchDirection() = 0;

  inline bool getStatusFlag() {
    return fStatusFlag;
  }

  virtual void reset();

  /** @brief See if the track representation has auxiliary information stored.
   * 
   * See if auxiliary information is stored in the track representation. See the
   * documentation of GFAbsTrackRep::getAuxInfo() for details.
   */
  virtual bool hasAuxInfo() {
    return false;
  }

  /** @brief Get auxiliary information from the track representation.
   *
   * AuxInfo is a mechanism which allows creators of track representations to
   * hand out any information they like (as long as it is compatible with a
   * TMatrixD). It should be used if setData requires additional 
   * information to update the representation, but it can also be used for
   * debugging information if needed. See also the documentation of 
   * GFAbsTrackRep::setData().
   */
  virtual const TMatrixD* getAuxInfo(const GFDetPlane&) {
    return NULL;
  }


  void resetXX0() {
    fXX0 = 0.;
  }

 private:
  void Abort(std::string method);


  /*----- Data members -----*/
 protected:
  //! Dimensionality of track representation
  unsigned int fDimension;

  //! The vector of track parameters
  TVectorD fState;

  //! The covariance matrix
  TMatrixDSym fCov;

  //! chiSqu of the track fit
  double       fChiSqu;
  double       fForwardChiSqu;
  double	   fNdf; // DAF needs non integer ndf

  //! status of track representation: 0 means everything's OK
  int fStatusFlag;

  //!state, cov and plane for first and last point in fit
  TVectorD    fFirstState;
  TMatrixDSym fFirstCov;

  TVectorD    fLastState;
  TMatrixDSym fLastCov;
  GFDetPlane  fFirstPlane;
  GFDetPlane  fLastPlane;

  // detector plane where the track parameters are given
  GFDetPlane  fRefPlane;

  double fXX0;


  ClassDef(GFAbsTrackRep,5)

};

#endif

/** @} */