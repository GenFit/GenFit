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

#ifndef GFDETPLANE_H
#define GFDETPLANE_H

#include "GFAbsFinitePlane.h"

#include "TObject.h"
#include "TVector3.h"

class TPolyMarker3D;
class TPolyLine3D;

/** @brief Detector plane genfit geometry class
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

class GFDetPlane : public TObject {
public:

  // Constructors/Destructors ---------
  GFDetPlane(GFAbsFinitePlane* finite=NULL);

  GFDetPlane(const TVector3& o,
             const TVector3& u,
             const TVector3& v,
             GFAbsFinitePlane* finite=NULL);

  GFDetPlane(const TVector3& o,
             const TVector3& n,
             GFAbsFinitePlane* finite=NULL);

  virtual ~GFDetPlane();

  GFDetPlane(const GFDetPlane&);

  GFDetPlane& operator=(const GFDetPlane&);

  // Accessors -----------------------
  const TVector3& getO() const {return fO;}
  const TVector3& getU() const {return fU;}
  const TVector3& getV() const {return fV;}

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
  void setFinitePlane(GFAbsFinitePlane* finite){fFinitePlane=finite;}

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
  TVector2 straightLineToPlane(const TVector3& point,const TVector3& dir) const;


  void Print(const Option_t* = "") const;

  //! for poor attempts of making an event display. There is a lot of room for improvements.
  void getGraphics(double mesh, double length, TPolyMarker3D **pl, TPolyLine3D **plLine,TPolyLine3D **u, TPolyLine3D **v, TPolyLine3D **n=NULL) const;

  //! this operator is called very often in Kalman filtering. It checks equality of planes
  //! by comparing the 9 double values that define them.
  friend bool operator== (const GFDetPlane& lhs, const GFDetPlane& rhs);
  //! returns NOT ==
  friend bool operator!= (const GFDetPlane& lhs, const GFDetPlane& rhs);

  //! absolute distance from a point to the plane
  double distance(const TVector3& point) const;
  double distance(double, double, double) const;


  //! intersect in the active area? C.f. GFAbsFinitePlane
  bool inActive(const TVector3& point, const TVector3& dir) const {
    if(fFinitePlane==NULL) return true;
    return this->inActive( this->straightLineToPlane(point,dir));
  }

  //! inActive methods refer to finite plane. C.f. GFAbsFinitePlane
  bool inActive(double u, double v) const{
    if(fFinitePlane==NULL) return true;
    return fFinitePlane->inActive(u,v);
  }

  //! inActive methods refer to finite plane. C.f. GFAbsFinitePlane
  bool inActive(const TVector2& v) const{
    return inActive(v.X(),v.Y());
  }
  
  bool isFinite() const {
    if (fFinitePlane != NULL) return true;
    return false;
  }

  // delete fFinitePlane and set O, U, V to default values
  void reset();


 private:

  // Private Methods -----------------
  void sane(); // ensures orthonormal coordinates

  // Private Data Members ------------
  // origin
  TVector3 fO;
  // Vectors spanning the plane
  TVector3 fU;
  TVector3 fV;

  GFAbsFinitePlane* fFinitePlane;

public:
  ClassDef(GFDetPlane,2)

};

#endif

/** @} */
