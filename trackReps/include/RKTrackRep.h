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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_RKTrackRep_h
#define genfit_RKTrackRep_h

#include "AbsTrackRep.h"
#include "StateOnPlane.h"
#include "StepLimits.h"
#include "Material.h"
#include "EigenMatrixTypedefs.h"

#include <algorithm>

namespace genfit {

/**
 * @brief Helper for RKTrackRep
 */
struct RKStep {
  MatStep matStep_; // material properties and stepsize
  Vector7 state7_; // 7D state vector
  StepLimits limits_;

  RKStep() : state7_(Vector7::Zero()) {}
};


/**
 * @brief Helper for RKTrackRep
 */
struct ExtrapStep {
  Matrix7x7 jac7_; // 5D jacobian of transport
  Matrix7x7Sym noise7_; // 5D noise matrix

  ExtrapStep() : jac7_(Matrix7x7::Zero()), noise7_(Matrix7x7Sym::Zero()) {}
};


/**
 * @brief AbsTrackRep with 5D track parameterization in plane coordinates: (q/p, u', v', u, v)
 *
 * q/p is charge over momentum.
 * u' and v' are direction tangents.
 * u and v are positions on a DetPlane.
 */
class RKTrackRep : public AbsTrackRep {
    friend class RKTrackRepTests_momMag_Test;
    friend class RKTrackRepTests_calcForwardJacobianAndNoise_Test;
    friend class RKTrackRepTests_getState7_Test;
    friend class RKTrackRepTests_getState5_Test;
    friend class RKTrackRepTests_calcJ_Mp_7x5_Test;
    friend class RKTrackRepTests_calcJ_pM_5x7_Test;

 public:

  RKTrackRep();
  RKTrackRep(int pdgCode, char propDir = 0);

  virtual ~RKTrackRep();

  virtual AbsTrackRep* clone() const {return new RKTrackRep(*this);}

  virtual double extrapolateToPlane(StateOnPlane& state,
      const SharedPlanePtr& plane,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  using AbsTrackRep::extrapolateToLine;

  virtual double extrapolateToLine(StateOnPlane& state,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const {
    return extrapToPoint(state, point, nullptr, stopAtBoundary, calcJacobianNoise);
  }

  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      const TMatrixDSym& G, // weight matrix (metric)
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const {
    return extrapToPoint(state, point, &G, stopAtBoundary, calcJacobianNoise);
  }

  virtual double extrapolateToCylinder(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  
  virtual double extrapolateToCone(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  virtual double extrapolateToSphere(StateOnPlane& state,
      double radius,
      const TVector3& point = TVector3(0.,0.,0.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  virtual double extrapolateBy(StateOnPlane& state,
      double step,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;


  unsigned int getDim() const {return 5;}

  virtual TVector3 getPos(const StateOnPlane& state) const;

  virtual TVector3 getMom(const StateOnPlane& state) const;
  virtual void getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const;

  virtual double getMomMag(const StateOnPlane& state) const;
  virtual double getMomVar(const MeasuredStateOnPlane& state) const;

  virtual TMatrixDSym get6DCov(const MeasuredStateOnPlane& state) const;
  virtual void getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const;
  virtual double getCharge(const StateOnPlane& state) const;
  virtual double getQop(const StateOnPlane& state) const {return state.getState()(0);}
  double getSpu(const StateOnPlane& state) const;
  double getTime(const StateOnPlane& state) const;

  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const;

  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const;

  std::vector<genfit::MatStep> getSteps() const;

  virtual double getRadiationLenght() const;

  virtual void setPosMom(StateOnPlane& state, const TVector3& pos, const TVector3& mom) const;
  virtual void setPosMom(StateOnPlane& state, const TVectorD& state6) const;
  virtual void setPosMomErr(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const;
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const;
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVectorD& state6, const TMatrixDSym& cov6x6) const;

  virtual void setChargeSign(StateOnPlane& state, double charge) const;
  virtual void setQop(StateOnPlane& state, double qop) const {state.getState()(0) = qop;}

  void setSpu(StateOnPlane& state, double spu) const;
  void setTime(StateOnPlane& state, double time) const;

  //! The actual Runge Kutta propagation
  /** propagate state7 with step S. Fills SA (Start directions derivatives dA/S).
   *  This is a single Runge-Kutta step.
   *  If jacobian is nullptr, only the state is propagated,
   *  otherwise also the 7x7 jacobian is calculated.
   *  If varField is false, the magnetic field will only be evaluated at the starting position.
   *  The return value is an estimation on how good the extrapolation is, and it is usually fine if it is > 1.
   *  It gives a suggestion how you must scale S so that the quality will be sufficient.
   */
  double RKPropagate(Vector7& state7,
                     Matrix7x7* jacobian,
                     Vector3& SA,
                     double S,
                     bool varField = true,
                     bool calcOnlyLastRowOfJ = false) const;
  virtual bool isSameType(const AbsTrackRep* other);
  virtual bool isSame(const AbsTrackRep* other);

 private:

  void initArrays() const;

  virtual double extrapToPoint(StateOnPlane& state,
      const TVector3& point,
      const TMatrixDSym* G = nullptr, // weight matrix (metric)
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  /***
   * Get 7d representation of a state on plane.
   * @param state - 5d state
   * @return 7d state
   */
  Vector7 getState7(const StateOnPlane& state) const;

  /***
   * Get 5d representation of a 7d state.
   * The 7d state must already like on plane of state
   * @param state - plane on which the state is projected
   * @param state7 - 7d state
   */
  void getState5(StateOnPlane& state, const Vector7& state7) const;

  Matrix5x7 calcJ_pM_5x7(const Vector7& state7, const DetPlane& plane) const;

  Matrix7x5 calcJ_Mp_7x5(const Vector7& state7, const DetPlane& plane) const;

  void calcForwardJacobianAndNoise(const Vector7& startState7, const DetPlane& startPlane,
                                   const Vector7& destState7, const DetPlane& destPlane) const;

  /***
   * Transform state6 covariance to state7 covariance.
   * ! plane and charge must already be set!
   * @param cov
   * @param state7
   * @param state
   */
  void transformM6P(const Matrix6x6Sym& cov, const Vector7& state7, MeasuredStateOnPlane& state) const;

  /***
   * Transform state5 covariance to state6 covariance.
   * @param state
   * @param out6x6
   */
  Matrix6x6Sym transformPM6(const MeasuredStateOnPlane& state) const;


  //! Propagates the particle through the magnetic field.
  /** If the propagation is successful and the plane is reached, the function returns true.
    * Propagated state and the jacobian of the extrapolation are written to state7 and jacobianT.
    * The jacobian is only calculated if jacobianT != nullptr.
    * In the main loop of the Runge Kutta algorithm, the estimateStep() is called
    * and may reduce the estimated stepsize so that a maximum momentum loss will not be exceeded,
    * and stop at material boundaries.
    * If this is the case, RKutta() will only propagate the reduced distance and then return. This is to ensure that
    * material effects, which are calculated after the propagation, are taken into account properly.
    */
  bool RKutta(const Vector4& SU,
              const DetPlane& plane,
              double charge,
              double mass,
              Vector7& state7,
              Matrix7x7* jacobianT,
              Vector7* J_MMT_unprojected_lastRow,
              double& coveredDistance, // signed
              double& flightTime,
              bool& checkJacProj,
              Matrix7x7Sym& noiseProjection,
              StepLimits& limits,
              bool onlyOneStep = false,
              bool calcOnlyLastRowOfJ = false) const;

  double estimateStep(const Vector7& state7,
                      const Vector4& SU,
                      const DetPlane& plane,
                      const Scalar& charge,
                      Scalar& relMomLoss,
                      StepLimits& limits) const;

  TVector3 pocaOnLine(const TVector3& linePoint,
                     const TVector3& lineDirection,
                     const TVector3& point) const;

  //! Handles propagation and material effects
  /** #extrapolateToPlane(), #extrapolateToPoint() and #extrapolateToLine() etc. call this function.
    * #Extrap() needs a plane as an argument, hence #extrapolateToPoint() and #extrapolateToLine() create virtual detector planes.
    * In this function, #RKutta() is called and the resulting points and point paths are filtered
    * so that the direction doesn't change and tiny steps are filtered out.
    * After the propagation the material effects are called via the MaterialEffects singleton.
    * #Extrap() will loop until the plane is reached, unless the propagation fails or the maximum number of
    * iterations is exceeded.
    */
  double Extrap(const DetPlane& startPlane, // plane where Extrap starts
                const DetPlane& destPlane, // plane where Extrap has to extrapolate to
                double charge,
                double mass,
                bool& isAtBoundary,
                Vector7& state7,
                double& flightTime,
                bool fillExtrapSteps,
                TMatrixDSym* cov = nullptr,
                bool onlyOneStep = false,
                bool stopAtBoundary = false,
                double maxStep = 1.E99) const;

  void checkCache(const StateOnPlane& state, const SharedPlanePtr* plane) const;

  mutable StateOnPlane lastStartState_; //! state where the last extrapolation has started
  mutable StateOnPlane lastEndState_; //! state where the last extrapolation has ended
  mutable std::vector<RKStep> RKSteps_; //! RungeKutta steps made in the last extrapolation
  mutable int RKStepsFXStart_; //!
  mutable int RKStepsFXStop_; //!
  mutable std::vector<ExtrapStep> ExtrapSteps_; //! steps made in Extrap during last extrapolation

  mutable Matrix5x5 fJacobian_; //!
  mutable Matrix5x5Sym fNoise_; //!

  mutable bool useCache_; //! use cached RKSteps_ for extrapolation
  mutable unsigned int cachePos_; //!

  // auxiliary variables and arrays
  // needed in Extrap()
  mutable StepLimits limits_; //!
  mutable Matrix7x7Sym noiseArray_; //! noise matrix of the last extrapolation
  mutable Matrix7x7Sym noiseProjection_; //!
  mutable Matrix7x7 J_MMT_; //!

 public:

  ClassDef(RKTrackRep, 2)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTrackRep_h
