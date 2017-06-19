/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

#ifndef genfit_AbsTrackRep_h
#define genfit_AbsTrackRep_h

#include "SharedPlanePtr.h"
//#include "MaterialInfo.h"
#include "Material.h"

#include <TVector3.h>
#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>


namespace genfit {

/**
 * @brief Simple struct containing MaterialProperties and stepsize in the material.
 */
struct MatStep {
  Material material_;
  double stepSize_;

  MatStep() {
    stepSize_ = 0;
  }

};

class StateOnPlane;
class MeasuredStateOnPlane;
class AbsMeasurement;

/**
 * @brief Abstract base class for a track representation
 *
 *  Provides functionality to extrapolate a StateOnPlane to another DetPlane,
 *  to the POCA to a line or a point, or a cylinder or sphere.
 *  Defines a set of parameters describing the track.
 *  StateOnPlane objects are always defined with a track parameterization of a specific AbsTrackRep.
 *  The AbsTrackRep provides functionality to translate from the internal representation of a state
 *  into cartesian position and momentum (and covariance) and vice versa.
 */
class AbsTrackRep : public TObject {

 public:

  AbsTrackRep();
  AbsTrackRep(int pdgCode, char propDir = 0);

  virtual ~AbsTrackRep() {;}

  //! Clone the trackRep.
  virtual AbsTrackRep* clone() const = 0;

  /**
   * @brief Extrapolates the state to plane, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToPlane(
      StateOnPlane& state,
      const genfit::SharedPlanePtr& plane,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Extrapolates the state to the POCA to a line, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToLine(StateOnPlane& state,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Resembles the interface of GFAbsTrackRep in old versions of genfit
   *
   * This interface to extrapolateToLine is intended to resemble the
   * interface of GFAbsTrackRep in old versions of genfit and is
   * implemented by default via the preceding function.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToLine(StateOnPlane& state,
      const TVector3& point1,
      const TVector3& point2,
      TVector3& poca,
      TVector3& dirInPoca,
      TVector3& poca_onwire,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const {
    TVector3 wireDir(point2 - point1);
    wireDir.Unit();
    double retval = this->extrapolateToLine(state, point1, wireDir, stopAtBoundary, calcJacobianNoise);
    poca = this->getPos(state);
    dirInPoca = this->getMom(state);
    dirInPoca.Unit();

    poca_onwire = point1 + wireDir*((poca - point1)*wireDir);
    
    return retval;
  }

  /**
   * @brief Extrapolates the state to the POCA to a point, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Extrapolates the state to the POCA to a point in the metric of G, and returns the extrapolation length
   *        and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      const TMatrixDSym& G, // weight matrix (metric)
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Extrapolates the state to the cylinder surface, and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToCylinder(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Extrapolates the state to the cone surface, and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToCone(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Extrapolates the state to the sphere surface, and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToSphere(StateOnPlane& state,
      double radius,
      const TVector3& point = TVector3(0.,0.,0.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  /**
   * @brief Extrapolates the state by step (cm) and returns the extrapolation length
   *       and, via reference, the extrapolated state.
   *
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   *
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateBy(StateOnPlane& state,
      double step,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  //! extrapolate to an AbsMeasurement
  double extrapolateToMeasurement(StateOnPlane& state,
      const AbsMeasurement* measurement,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  //! Get the dimension of the state vector used by the track representation.
  virtual unsigned int getDim() const = 0;

  //! Get the cartesian position of a state.
  virtual TVector3 getPos(const StateOnPlane& state) const = 0;

  //! Get the cartesian momentum vector of a state.
  virtual TVector3 getMom(const StateOnPlane& state) const = 0;

  //! Get the direction vector of a state.
  TVector3 getDir(const StateOnPlane& state) const {return getMom(state).Unit();}

  //! Get cartesian position and momentum vector of a state.
  virtual void getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const = 0;

  //! Get cartesian position and direction vector of a state.
  void getPosDir(const StateOnPlane& state, TVector3& pos, TVector3& dir) const {getPosMom(state, pos, dir); dir.SetMag(1.);}

  //! Get the 6D state vector (x, y, z, p_x, p_y, p_z).
  virtual TVectorD get6DState(const StateOnPlane& state) const;

  //! Get the 6D covariance.
  virtual TMatrixDSym get6DCov(const MeasuredStateOnPlane& state) const = 0;

  //! Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance.
  virtual void getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const = 0;

  //! Translates MeasuredStateOnPlane into 6D state vector (x, y, z, p_x, p_y, p_z) and 6x6 covariance.
  virtual void get6DStateCov(const MeasuredStateOnPlane& state, TVectorD& stateVec, TMatrixDSym& cov) const;

  //! get the magnitude of the momentum in GeV.
  virtual double getMomMag(const StateOnPlane& state) const = 0;
  //! get the variance of the absolute value of the momentum .
  virtual double getMomVar(const MeasuredStateOnPlane& state) const = 0;

  //! Get the pdg code.
  int getPDG() const {return pdgCode_;}

  //! Get the charge of the particle of the pdg code
  double getPDGCharge() const;

  /**
   * @brief Get the (fitted) charge of a state.
   * This is not always equal the pdg charge (e.g. if the charge sign was flipped during the fit).
   */
  virtual double getCharge(const StateOnPlane& state) const = 0;
  //! Get charge over momentum.
  virtual double getQop(const StateOnPlane& state) const = 0;
  //! Get tha particle mass in GeV/c^2
  double getMass(const StateOnPlane& state) const;

  //! Get propagation direction. (-1, 0, 1) -> (backward, auto, forward).
  char getPropDir() const {return propDir_;}

  //! Get the jacobian and noise matrix of the last extrapolation.
  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const = 0;

  //! Get the jacobian and noise matrix of the last extrapolation if it would have been done in opposite direction.
  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const = 0;

  //! Get stepsizes and material properties of crossed materials of the last extrapolation.
  virtual std::vector<genfit::MatStep> getSteps() const = 0;

  //! Get the accumulated X/X0 (path / radiation length) of the material crossed in the last extrapolation.
  virtual double getRadiationLenght() const = 0;

  //! Get the time corresponding to the StateOnPlane.  Extrapolation
  // should keep this up to date with the time of flight.
  virtual double getTime(const StateOnPlane&) const = 0;

  /**
   * @brief Calculate Jacobian of transportation numerically.
   * Slow but accurate. Can be used to validate (semi)analytic calculations.
   */
  void calcJacobianNumerically(const genfit::StateOnPlane& origState,
                                   const genfit::SharedPlanePtr destPlane,
                                   TMatrixD& jacobian) const;

  //! try to multiply pdg code with -1. (Switch from particle to anti-particle and vice versa).
  bool switchPDGSign();

  //! Set position and momentum of state.
  virtual void setPosMom(StateOnPlane& state, const TVector3& pos, const TVector3& mom) const = 0;
  //! Set position and momentum of state.
  virtual void setPosMom(StateOnPlane& state, const TVectorD& state6) const = 0;
  //! Set position and momentum and error of state.
  virtual void setPosMomErr(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const = 0;
  //! Set position, momentum and covariance of state.
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const = 0;
  //! Set position, momentum and covariance of state.
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVectorD& state6, const TMatrixDSym& cov6x6) const = 0;

  //! Set the sign of the charge according to charge.
  virtual void setChargeSign(StateOnPlane& state, double charge) const = 0;
  //! Set charge/momentum.
  virtual void setQop(StateOnPlane& state, double qop) const = 0;
  //! Set time at which the state was defined
  virtual void setTime(StateOnPlane& state, double time) const = 0;

  //! Set propagation direction. (-1, 0, 1) -> (backward, auto, forward).
  void setPropDir(int dir) {
    if (dir>0) propDir_ = 1;
    else if (dir<0) propDir_ = -1;
    else propDir_ = 0;
  };

  //! Switch propagation direction. Has no effect if propDir_ is set to 0.
  void switchPropDir(){propDir_ = -1*propDir_;}

  //! check if other is of same type (e.g. RKTrackRep).
  virtual bool isSameType(const AbsTrackRep* other) = 0;

  //! check if other is of same type (e.g. RKTrackRep) and has same pdg code.
  virtual bool isSame(const AbsTrackRep* other) = 0;

  virtual void setDebugLvl(unsigned int lvl = 1) {debugLvl_ = lvl;}

  virtual void Print(const Option_t* = "") const;

 protected:

  //! protect from calling copy c'tor from outside the class. Use #clone() if you want a copy!
  AbsTrackRep(const AbsTrackRep&);
  //! protect from calling assignment operator from outside the class. Use #clone() instead!
  AbsTrackRep& operator=(const AbsTrackRep&);


  //! Particle code
  int pdgCode_;
  //! propagation direction (-1, 0, 1) -> (backward, auto, forward)
  char propDir_;

  unsigned int debugLvl_;

 public:
  ClassDef(AbsTrackRep,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsTrackRep_h
