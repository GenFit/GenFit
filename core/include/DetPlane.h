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
// Description:
//      Detector plane - a geometric object
/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup genfit
 * @{
 */

#ifndef genfit_DetPlane_h
#define genfit_DetPlane_h

#include "AbsFinitePlane.h"

#include <TObject.h>
#include <TVector3.h>

#include <memory>


namespace genfit {

/** @brief Detector plane.
 *
 * A detector plane is the principle object to define coordinate systems for
 * track fitting in genfit. Since a particle trajectory is a
 * one-dimensional object (regardless of any specific parameterization)
 * positions with respect to the track are always measured in a plane.
 *
 * Which plane is chosen depends on the type of detector. Fixed plane
 * detectors have their detector plane defined by their mechanical setup. While
 * wire chambers or time projection chambers might want to define a detector
 * plane more flexibly.
 *
 * This class parameterizes a plane in terms of an origin vector o
 * and two plane-spanning directions u and v.
 */
class DetPlane : public TObject {

 public:


  // Constructors/Destructors ---------
  DetPlane(AbsFinitePlane* finite = nullptr);

  DetPlane(const TVector3& o,
             const TVector3& u,
             const TVector3& v,
             AbsFinitePlane* finite = nullptr);

  DetPlane(const TVector3& o,
             const TVector3& n,
             AbsFinitePlane* finite = nullptr);

  virtual ~DetPlane();

  DetPlane(const DetPlane&);
  DetPlane& operator=(DetPlane);
  void swap(DetPlane& other); // nothrow

  // Accessors -----------------------
  const TVector3& getO() const {return o_;}
  const TVector3& getU() const {return u_;}
  const TVector3& getV() const {return v_;}

  // Modifiers -----------------------
  void set(const TVector3& o,
           const TVector3& u,
           const TVector3& v);
  void setO(const TVector3& o);
  void setO(double, double, double);
  void setU(const TVector3& u);
  void setU(double, double, double);
  void setV(const TVector3& v);
  void setV(double, double, double);
  void setUV(const TVector3& u, const TVector3& v);
  void setON(const TVector3& o, const TVector3& n);

  //! Optionally, set the finite plane definition. This is most important for
  //! avoiding fake intersection points in fitting of curlers. This should
  //! be implemented for silicon detectors most importantly.
  void setFinitePlane(AbsFinitePlane* finite){finitePlane_.reset(finite);}

  // Operations ----------------------
  TVector3 getNormal() const;
  void setNormal(const TVector3& n);
  void setNormal(double, double, double);
  void setNormal(const double& theta, const double& phi);

  //! projecting a direction onto the plane:
  TVector2 project(const TVector3& x) const;

  //! transform from Lab system into plane
  TVector2 LabToPlane(const TVector3& x) const;

  //! transform from plane coordinates to lab system
  TVector3 toLab(const TVector2& x) const;

  // get vector from point to plane (normal)
  TVector3 dist(const TVector3& point) const;

  //! gives u,v coordinates of the intersection point of a straight line with plane
  TVector2 straightLineToPlane(const TVector3& point, const TVector3& dir) const;

  //! gives u,v coordinates of the intersection point of a straight line with plane
  void straightLineToPlane(const double& posX, const double& posY, const double& posZ,
                           const double& dirX, const double& dirY, const double& dirZ,
                           double& u, double& v) const;

  void Print(const Option_t* = "") const;

  //! Checks equality of planes by comparing the 9 double values that define them.
  friend bool operator== (const DetPlane& lhs, const DetPlane& rhs);
  //! returns NOT ==
  friend bool operator!= (const DetPlane& lhs, const DetPlane& rhs);

  //! absolute distance from a point to the plane
  double distance(const TVector3& point) const;
  double distance(double, double, double) const;


  //! intersect in the active area? C.f. AbsFinitePlane
  bool isInActive(const TVector3& point, const TVector3& dir) const {
    if(finitePlane_.get() == nullptr) return true;
    return this->isInActive( this->straightLineToPlane(point,dir));
  }

  //! intersect in the active area? C.f. AbsFinitePlane
  bool isInActive(const double& posX, const double& posY, const double& posZ,
                  const double& dirX, const double& dirY, const double& dirZ) const {
    if(finitePlane_.get() == nullptr) return true;
    double u, v;
    this->straightLineToPlane(posX, posY, posZ, dirX, dirY, dirZ, u, v);
    return this->isInActive(u, v);
  }

  //! isInActive methods refer to finite plane. C.f. AbsFinitePlane
  bool isInActive(double u, double v) const{
    if(finitePlane_.get() == nullptr) return true;
    return finitePlane_->isInActive(u,v);
  }

  //! isInActive methods refer to finite plane. C.f. AbsFinitePlane
  bool isInActive(const TVector2& v) const{
    return isInActive(v.X(),v.Y());
  }

  bool isFinite() const {
    return (finitePlane_.get() != nullptr);
  }

  //! rotate u and v around normal. Angle is in rad. More for debugging than for actual use.
  void rotate(double angle);

  //! delete finitePlane_ and set O, U, V to default values
  void reset();

 private:
  // Private Methods -----------------
  //! ensures orthonormal coordinates
  void sane();

  TVector3 o_;
  TVector3 u_;
  TVector3 v_;

  std::unique_ptr<AbsFinitePlane> finitePlane_; // Ownership

 public:
  ClassDef(DetPlane,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_DetPlane_h
