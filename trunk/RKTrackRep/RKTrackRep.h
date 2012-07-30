/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

/*
 */

/** @addtogroup RKTrackRep
 * @{
 */


#ifndef RKTRACKREP_H
#define RKTRACKREP_H


#include "GFAbsTrackRep.h"
#include "GFDetPlane.h"
#include "GFTrackCand.h"
#include <TMatrixD.h>

/** @brief Track Representation module based on a Runge-Kutta algorithm including a full material model
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, author)
 *
 * The Runge Kutta implementation stems from GEANT3 originally (R. Brun et al.).
 * Porting to %C goes back to Igor Gavrilenko @ CERN.
 * The code was taken from the Phast analysis package of the COMPASS experiment
 * (Sergei Gerrassimov @ CERN).
 *
 * The state is 5D: (q/p, u', v', u,v)
 */


class RKTrackRep : public GFAbsTrackRep {

 public:

  // Constructors/Destructors ---------
  RKTrackRep();
  RKTrackRep(const TVector3& pos,
       const TVector3& mom,
       const TVector3& poserr,
       const TVector3& momerr,
       const int& PDGCode);

  RKTrackRep(const TVector3& pos,
       const TVector3& mom,
       const int& PDGCode);

  RKTrackRep(const GFTrackCand* aGFTrackCandPtr);

  virtual ~RKTrackRep();


  virtual GFAbsTrackRep* clone() const {return new RKTrackRep(*this);}
  virtual GFAbsTrackRep* prototype()const{return new RKTrackRep();}

  //! returns the tracklength spanned in this extrapolation
  /** The covariance matrix is transformed from the plane coordinate system to the master reference system (for the propagation) and, after propagation, back to the plane coordinate system.\n
    * Also the parameter spu (which is +1 or -1 and indicates the direction of the particle) is calculated and stored in #fCacheSpu. The plane is stored in #fCachePlane.
    * \n
    *
    * Master reference system (MARS):
    * \f{eqnarray*}x & = & O_{x}+uU_{x}+vV_{x}\\y & = & O_{y}+uU_{y}+vV_{y}\\z & = & O_{z}+uU_{z}+vV_{z}\\a_{x} & = & \frac{\mbox{spu}}{\widetilde{p}}\left(N_{x}+u\prime U_{x}+v\prime V_{x}\right)\\a_{y} & = & \frac{\mbox{spu}}{\widetilde{p}}\left(N_{y}+u\prime U_{y}+v\prime V_{y}\right)\\a_{z} & = & \frac{\mbox{spu}}{\widetilde{p}}\left(N_{z}+u\prime U_{z}+v\prime V_{z}\right)\\\frac{q}{p} & = & \frac{q}{p}\f}
    * Plane coordinate system:
    * \f{eqnarray*}u & = & \left(x-O_{x}\right)U_{x}+\left(y-O_{y}\right)U_{y}+\left(z-O_{z}\right)U_{z}\\v & = & \left(x-O_{x}\right)V_{x}+\left(y-O_{y}\right)V_{y}+\left(z-O_{z}\right)V_{z}\\u\prime & = & \frac{a_{x}U_{x}+a_{y}U_{y}+a_{z}U_{z}}{a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}}\\v\prime & = & \frac{a_{x}V_{x}+a_{y}V_{y}+a_{z}V_{z}}{a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}}\\\frac{q}{p} & = & \frac{q}{p}\\\mbox{spu} & = & \frac{a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}}{|a_{x}^{2}N_{x}^{2}+a_{y}^{2}N_{y}^{2}+a_{z}^{2}N_{z}^{2}|}=\pm1\f}
    *
    * Jacobians:\n
    * \f$J_{p,M}=\left(\begin{array}{ccccc}\frac{\partial x}{\partial u} & \frac{\partial x}{\partial v} & \frac{\partial x}{\partial u\prime} & \frac{\partial x}{\partial v\prime} & \frac{\partial x}{\partial\frac{q}{p}}\\\frac{\partial y}{\partial u} & \frac{\partial y}{\partial v} & \frac{\partial y}{\partial u\prime} & \frac{\partial y}{\partial v\prime} & \frac{\partial y}{\partial\frac{q}{p}}\\\frac{\partial z}{\partial u} & \frac{\partial z}{\partial v} & \frac{\partial z}{\partial u\prime} & \frac{\partial z}{\partial v\prime} & \frac{\partial z}{\partial\frac{q}{p}}\\\frac{\partial a_{x}}{\partial u} & \frac{\partial a_{x}}{\partial v} & \frac{\partial a_{x}}{\partial u\prime} & \frac{\partial a_{x}}{\partial v\prime} & \frac{\partial a_{x}}{\partial\frac{q}{p}}\\\frac{\partial a_{y}}{\partial u} & \frac{\partial a_{y}}{\partial v} & \frac{\partial a_{y}}{\partial u\prime} & \frac{\partial a_{y}}{\partial v\prime} & \frac{\partial a_{y}}{\partial\frac{q}{p}}\\\frac{\partial a_{z}}{\partial u} & \frac{\partial a_{z}}{\partial v} & \frac{\partial a_{z}}{\partial u\prime} & \frac{\partial a_{z}}{\partial v\prime} & \frac{\partial a_{z}}{\partial\frac{q}{p}}\\\frac{\partial\frac{q}{p}}{\partial u} & \frac{\partial\frac{q}{p}}{\partial v} & \frac{\partial\frac{q}{p}}{\partial u\prime} & \frac{\partial\frac{q}{p}}{\partial v\prime} & \frac{\partial\frac{q}{p}}{\partial\frac{q}{p}}\end{array}\right)\f$
    *
    * \f$J_{p,M}=\left(\begin{array}{cccccc}U_{x} & V_{x} & 0 & 0 & 0\\U_{y} & V_{y} & 0 & 0 & 0\\U_{z} & V_{z} & 0 & 0 & 0\\0 & 0 & \left\{ \frac{\textrm{spu}}{\widetilde{p}}U_{x}-\left[\frac{\textrm{spu}}{\widetilde{p}^{3}}\widetilde{p}_{x}\left(U_{x}\widetilde{p}_{x}+U_{y}\widetilde{p}_{y}+U_{z}\widetilde{p}_{z}\right)\right]\right\}  & \left\{ \frac{\textrm{spu}}{\widetilde{p}}V_{x}-\left[\frac{\textrm{spu}}{\widetilde{p}^{3}}\widetilde{p}_{x}\left(V_{x}\widetilde{p}_{x}+V_{y}\widetilde{p}_{y}+V_{z}\widetilde{p}_{z}\right)\right]\right\}  & 0\\0 & 0 & \left\{ \frac{\textrm{spu}}{\widetilde{p}}U_{y}-\left[\frac{\textrm{spu}}{\widetilde{p}^{3}}\widetilde{p}_{y}\left(U_{x}\widetilde{p}_{x}+U_{y}\widetilde{p}_{y}+U_{z}\widetilde{p}_{z}\right)\right]\right\}  & \left\{ \frac{\textrm{spu}}{\widetilde{p}}V_{y}-\left[\frac{\textrm{spu}}{\widetilde{p}^{3}}\widetilde{p}_{y}\left(V_{x}\widetilde{p}_{x}+V_{y}\widetilde{p}_{y}+V_{z}\widetilde{p}_{z}\right)\right]\right\}  & 0\\0 & 0 & \left\{ \frac{\textrm{spu}}{\widetilde{p}}U_{z}-\left[\frac{\textrm{spu}}{\widetilde{p}^{3}}\widetilde{p}_{z}\left(U_{x}\widetilde{p}_{x}+U_{y}\widetilde{p}_{y}+U_{z}\widetilde{p}_{z}\right)\right]\right\}  & \left\{ \frac{\textrm{spu}}{\widetilde{p}}V_{z}-\left[\frac{\textrm{spu}}{\widetilde{p}^{3}}\widetilde{p}_{z}\left(V_{x}\widetilde{p}_{x}+V_{y}\widetilde{p}_{y}+V_{z}\widetilde{p}_{z}\right)\right]\right\}  & 0\\0 & 0 & 0 & 0 & 1\end{array}\right)\f$
    * 
    * with
    * \f{eqnarray*}\widetilde{p} & = & \sqrt{\widetilde{p}_{x}^{2}+\widetilde{p}_{y}^{2}+\widetilde{p}_{z}^{2}}\\\widetilde{p}_{x} & = & \textrm{spu}\left(N_{x}+u\prime U_{x}+v\prime V_{x}\right)\\\widetilde{p}_{y} & = & \textrm{spu}\left(N_{y}+u\prime U_{y}+v\prime V_{y}\right)\\\widetilde{p}_{z} & = & \textrm{spu}\left(N_{z}+u\prime U_{z}+v\prime V_{z}\right)\f}
    *
    * \f$J_{M,p}=\left(\begin{array}{ccccccc}\frac{\partial u}{\partial x} & \frac{\partial u}{\partial y} & \frac{\partial u}{\partial z} & \frac{\partial u}{\partial a_{x}} & \frac{\partial u}{\partial a_{y}} & \frac{\partial u}{\partial a_{z}} & \frac{\partial u}{\partial\frac{q}{p}}\\\frac{\partial v}{\partial x} & \frac{\partial v}{\partial y} & \frac{\partial v}{\partial z} & \frac{\partial v}{\partial a_{x}} & \frac{\partial v}{\partial a_{y}} & \frac{\partial v}{\partial a_{z}} & \frac{\partial v}{\partial\frac{q}{p}}\\\frac{\partial u\prime}{\partial x} & \frac{\partial u\prime}{\partial y} & \frac{\partial u\prime}{\partial z} & \frac{\partial u\prime}{\partial a_{x}} & \frac{\partial u\prime}{\partial a_{y}} & \frac{\partial u\prime}{\partial a_{z}} & \frac{\partial u\prime}{\partial\frac{q}{p}}\\\frac{\partial v\prime}{\partial x} & \frac{\partial v\prime}{\partial y} & \frac{\partial v\prime}{\partial z} & \frac{\partial v\prime}{\partial a_{x}} & \frac{\partial v\prime}{\partial a_{y}} & \frac{\partial v\prime}{\partial a_{z}} & \frac{\partial v\prime}{\partial\frac{q}{p}}\\ \frac{\partial\frac{q}{p}}{\partial x} & \frac{\partial\frac{q}{p}}{\partial y} & \frac{\partial\frac{q}{p}}{\partial z} & \frac{\partial\frac{q}{p}}{\partial a_{x}} & \frac{\partial\frac{q}{p}}{\partial a_{y}} & \frac{\partial\frac{q}{p}}{\partial a_{z}} & \frac{\partial\frac{q}{p}}{\partial\frac{q}{p}}\\ \\\end{array}\right)\f$ 
    *
    * \f$J_{M,p}=\left(\begin{array}{ccccccc} U_{x} & U_{y} & U_{z} & 0 & 0 & 0 & 0\\ V_{x} & V_{y} & V_{z} & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & \frac{U_{x}\left(a_{y}N_{y}+a_{z}N_{z}\right)-N_{x}\left(a_{y}U_{y}+a_{z}U_{z}\right)}{\left(a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}\right)^{2}} & \frac{U_{y}\left(a_{x}N_{x}+a_{z}N_{z}\right)-N_{y}\left(a_{x}U_{x}+a_{z}U_{z}\right)}{\left(a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}\right)^{2}} & \frac{U_{z}\left(a_{x}N_{x}+a_{y}N_{y}\right)-N_{z}\left(a_{x}U_{x}+a_{y}U_{y}\right)}{\left(a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}\right)^{2}} & 0\\ 0 & 0 & 0 & \frac{V_{x}\left(a_{y}N_{y}+a_{z}N_{z}\right)-N_{x}\left(a_{y}V_{y}+a_{z}V_{z}\right)}{\left(a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}\right)^{2}} & \frac{V_{y}\left(a_{x}N_{x}+a_{z}N_{z}\right)-N_{y}\left(a_{x}V_{x}+a_{z}V_{z}\right)}{\left(a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}\right)^{2}} & \frac{V_{z}\left(a_{x}N_{x}+a_{y}N_{y}\right)-N_{z}\left(a_{x}V_{x}+a_{y}V_{y}\right)}{\left(a_{x}N_{x}+a_{y}N_{y}+a_{z}N_{z}\right)^{2}} & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \\\end{array}\right)\f$
    *
    */
  double extrapolate(const GFDetPlane&, 
         TMatrixD& statePred,
         TMatrixD& covPred);

  //! returns the tracklength spanned in this extrapolation
  double extrapolate(const GFDetPlane&, 
         TMatrixD& statePred);

  //! This method is to extrapolate the track to point of closest approach to a point in space
  void extrapolateToPoint(const TVector3& pos,
         TVector3& poca,
         TVector3& dirInPoca);
  
  //! This method extrapolates to the point of closest approach to a line
  void extrapolateToLine(const TVector3& point1,
         const TVector3& point2,
         TVector3& poca,
         TVector3& dirInPoca,
         TVector3& poca_onwire);
  
  //! make step of h cm along the track, returns the tracklength spanned in this extrapolation
  /** Also returns the position and direction by reference.
  * It does NOT alter the state of the trackrep and starts extrapolating from #fRefPlane.
  */
  double stepalong(double h,
         TVector3& point,
         TVector3& dir);

  //! Returns position of the track in the plane  
  /** If #GFDetPlane equals the reference plane #fRefPlane, returns current position; otherwise it extrapolates 
    * the track to the plane and returns the position.
    */
  TVector3 getPos(const GFDetPlane&);
  //! Returns momentum of the track in the plane  
  /** If #GFDetPlane equals the reference plane #fRefPlane, returns current momentum; otherwise it extrapolates 
    * the track to the plane and returns the momentum.
    */
  TVector3 getMom(const GFDetPlane&);
  //! Gets position and momentum in the plane 
  /** If #GFDetPlane equals the reference plane #fRefPlane, it gets current position and momentum; otherwise it extrapolates 
    * the track to the plane and gets the position and momentum.
    */
  void getPosMom(const GFDetPlane&,TVector3& pos,TVector3& mom);

  void getPosMomCov(const GFDetPlane& pl,
                    TVector3& pos, TVector3& mom,
                    TMatrixD& cov);

  //! Returns charge
  double getCharge()const {return fCharge;}

  int getPDG() {return fPdg;};

  //! Set propagation direction. (-1,0,1) -> (backward prop,decide myself,forward)
  void setPropDir(int dir);

  //! Switch propagation direction. Has no effect if propdir is set to 0.
  void switchDirection(){fDirection = -1*fDirection;}

  //! Set PDG particle code
  void setPDG(int);

  //! Sets state, plane and (optionally) covariance
  /** This function also sets the parameter #fSpu to the value stored in #fCacheSpu. Therefore it has to be ensured that
    * the plane #pl is the same as the plane of the last extrapolation (i.e. #fCachePlane), where #fCacheSpu was calculated.
    * Hence, if the argument #pl is not equal to #fCachePlane, an error message is shown an an exception is thrown.
    */
  void setData(const TMatrixD& st,
               const GFDetPlane& pl,
               const TMatrixD* cov=NULL,
               const TMatrixD* aux=NULL);

  //! Sets state, plane and covariance from position, momentum and 6x6 covariance
  /** Also sets the reference plane at position
    */
  void setPosMomCov(const TVector3& pos,
                    const TVector3& mom,
                    const TMatrixD& cov);

  void disableMaterialEffects(bool opt = true){fNoMaterial = opt;}

  const TMatrixD* getAuxInfo(const GFDetPlane& pl);
  
  bool hasAuxInfo() { return true; }



 private:

  void calcStateCov(const TVector3& pos,
                    const TVector3& mom,
                    const TVector3& poserr,
                    const TVector3& momerr);

  void calcState(const TVector3& pos,
                 const TVector3& mom);

  TMatrixD getState7() const;
  TMatrixD getState7(const TMatrixD& state5, const GFDetPlane& pl, const double& spu) const;
  TMatrixD getState5(const TMatrixD& state7, const GFDetPlane& pl, double& spu) const;

  void transformPM(const TMatrixD& in5x5,
                   TMatrixD& out,
                   const GFDetPlane& pl,
                   const TMatrixD& state5,
                   const double& spu,
                   TMatrixD* Jac = NULL) const;

  void transformMP(const TMatrixD& in,
                   TMatrixD& out5x5,
                   const GFDetPlane& pl,
                   const TMatrixD& state7,
                   TMatrixD* Jac = NULL) const;

  RKTrackRep& operator=(const RKTrackRep* rhs){return *this;};

  //! Propagates the particle through the magnetic field.
  /** If the propagation is successfull and the plane is reached, the function returns true.
    * The argument P has to contain the state (#P[0] - #P[6]) and a unity matrix (#P[7] - #P[55]) 
    * with the last column multiplied wit q/p (hence #P[55] is not 1 but q/p).  
    * Propagated state and the jacobian (with the last column multiplied wit q/p) of the extrapolation are written to #P.
    * In the main loop of the Runge Kutta algorithm, the steppers in  #fEffect are called
    * and may reduce the estimated stepsize so that a maximum momentum loss will not be exceeded.
    * If this is the case, RKutta() will only propagate the reduced distance and then return. This is to ensure that 
    * material effects, which are calculated after the propagation, are taken into account properly.
    * 
    */
  bool RKutta (const GFDetPlane& plane,
               double* P,
               double& coveredDistance,
               std::vector<TVector3>& points,
               std::vector<double>& pointLengths,
               bool calcCov=true) const;

  double estimateStep(const TVector3& pos,
                      const TVector3& dir,
                      const double* SU,
                      const GFDetPlane& plane,
                      const double& mom,
                      double& relMomLoss,
                      double& deltaAngle,
                      bool& stopBecauseOfMaterial) const;

  TVector3 poca2Line(const TVector3& extr1,
                     const TVector3& extr2,
                     const TVector3& point) const;
    
  //! Handles propagation and material effects
  /** extrapolate(), extrapolateToPoint() and extrapolateToLine() call this function.
    * Extrap() needs a plane as an argument, hence extrapolateToPoint() and extrapolateToLine() create virtual detector planes.
    * In this function, RKutta() is called and the resulting points and point paths are filtered 
    * so that the direction doesn't change and tiny steps are filtered out. After the propagation the material effects are called via the GFMaterialEffects singleton.
    * Extrap() will loop until the plane is reached, unless the propagation fails or the maximum number of 
    * iterations is exceeded.
    * fXX0 is also updated here.
    */
  double Extrap(const GFDetPlane& plane,
                TMatrixD* state,
                TMatrixD* cov=NULL);
  
  //RKTrackRep(const RKTrackRep& rhs){};
  
  // data members
  
  int fDirection;   // (-1,0,1) -> (backward prop,decide myself,forward)
  bool fNoMaterial; // don't calculate material effects if true
    
  //! PDG particle code
  int fPdg;
  //! Mass (in GeV)
  double fMass;
  //! Charge
  double fCharge;

  GFDetPlane fCachePlane;
  double fCacheSpu;
  double fSpu;
  TMatrixD fAuxInfo;

 public:
  ClassDef(RKTrackRep,5)

};

#endif

/** @} */
