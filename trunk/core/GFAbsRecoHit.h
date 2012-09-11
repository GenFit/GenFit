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
/** @addtogroup genfit */
/* @{ */



#ifndef GFABSRECOHIT_H
#define GFABSRECOHIT_H

#include<iostream>

#include "TMatrixT.h"
#include "TObject.h"

#include "GFAbsTrackRep.h"
#include "GFDetPlane.h"

/** @brief Base Class for representing a Hit in GENFIT
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * A hit is defined as a single measurement of a detector. Each detector can 
 * define it's own hit representation (geometry) by inherting from GFAbsRecoHit. 
 * We call such a child object a "RecoHit" to make clear that
 * it inherits from GFAbsRecoHit.
 * All detector specific information is handled inside the RecoHit objects.
 * The idea behind this is that the fitting algorithms can work on different
 * detectors simultanously. 
 *
 * GFAbsRecoHit defines the basic interface that is used by all genfit algorithms
 * to access hit-measurement information. It provides:
 *  - the matrices to store hit coordinates and covariances
 *  - defines methods to access these data
 *  - the interface to access a hit's detector plane object
 * 
 * All hits have to inherit from this base class. 
 * Inheritance can be direct or through template class 
 * RecoHitIfc<GeometryPolicy>.
 * These interface classes (defined with a specific policy)
 * provide additional functionality for specific hit geometries, 
 * such as space points, wires, etc. For details look 
 * at the RecoHitIfc documentation.
 * 
 * A simple example is given in VirtSpacePointRecoHit
 *
 * Background information: The main feature here is
 * that coordinates and covariances are available as general 
 * TMatrixT<double> objects. By using these general algebraic 
 * objects it is possible to abstract from the concrete measurement and 
 * thus define a unified framework to handle the data. This is a key ingredient
 * for the separation of data and algorithms which in turn enables us 
 * to elegantly combine information from different detectors.
 */

class GFAbsRecoHit : public TObject{
protected:
  /// Vector of raw coordinates of hit
  TMatrixT<double> fHitCoord;

  /// Covariance of raw hit coordinates
  TMatrixT<double> fHitCov;

  /// Sorting parameter used in GFTrack::sortHits()
  double fRho;

 private:
  int fNparHit;

public:
  virtual ~GFAbsRecoHit();

  /** @brief Constructor specifying dimension of coordinate vector
   *
   * One hit is generally represented by a vector of coordinates. 
   * @param NparHit specifies the dimension of this vector. 
   * (e.g. 3 for a spacepoint, 2 for a pixel, 1 for a strip)
   */
  GFAbsRecoHit(int NparHit);

  /** @brief Default constructor needed for compatibility with ROOT
   */
  GFAbsRecoHit();

  /** @brief Get transformation matrix. Transformation between hit 
   * coordinates and track representation coordinates.
   *
   * This is a virtual abstract method which has to be implemented in the child
   * classes.
   *
   * In general there is a linear transformation between the coordinate system
   * of the hit (which is defined by the detector plane) and the coordinates of
   * the track representation in that plane. In the most simple case the track
   * representation has 5 parameters (space + momentum) 
   * while a hit usually has less (one to three space coordinates). 
   *
   * The transformation matrix is then simply projecting out the 
   * space-components of the track representation.
   *
   * Its dimensions are NxM. Where N is the number of dimensions of the 
   * hit in the detector plane (usually 2 or 1) and M is the dimension of the
   * track representation.
   *
   * In this method a hit has to define with which track representations it can
   * work together. It should be the only point where this explicit 
   * coordination is necessary. 
   *
   * For example code see implementing classes below:
   */
  virtual TMatrixT<double> getHMatrix(const GFAbsTrackRep* stateVector)=0;

  /** @brief get measurement vector and hit covariance
   *
   */
  virtual void getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TMatrixT<double>& statePred,const TMatrixT<double>& covPred,TMatrixT<double>& m, TMatrixT<double>& V) = 0;
  
  /** @brief Get raw hit covariances. 
   *
   */
  TMatrixT<double> getRawHitCov() const {return fHitCov;}

  /** @brief Get raw hit coordinates. 
   *
   */
  TMatrixT<double> getRawHitCoord() const {return fHitCoord;}
  
  
  /** @brief Get detector plane for a given track representation.
   *
   * Virtual abstract method has to be implemented by inherting classes.
   *
   * In general the detector plane can depend both on the detector/hit geometry
   * as well as the track geometry. E.g. for a space point one usually chooses 
   * a plane that is perpendicular to the current track, since in that case no 
   * other plane is predefined.
   *
   * There are several implementations for this method in the HitPolicies. 
   * In the most simple case (a planar detector) the method just returns a
   * fixed (detector module specific) plane. This behaviour for example is 
   * implemented in PlanarHitPolicy.
   */
  virtual const GFDetPlane& getDetPlane(GFAbsTrackRep*) = 0; 
   
  /**
   * Get sorting parameter rho.
   */
  double getRho() {return fRho;}

  /**
   * Set sorting parameter rho.
   */
  void setRho(double rho) {fRho = rho;}
  
  /** @brief Get clone of this object.
   *
   * Virtual abstract method. Has to be implemented by inherting classes.
   * Creates a deep copy of this object. 
   * Ownership is trandsferred to the caller!
   */
  virtual GFAbsRecoHit* clone() = 0;

  /** @brief Print raw hit coordinates.
   */
  virtual void Print(const Option_t* option = "") const {fHitCoord.Print(option);}

  virtual const std::string& getPolicyName();

  int getNparHit(){return fNparHit;}

 public:
  ClassDef(GFAbsRecoHit,4)

};
  

#endif //FITTER_ABSHIT_H

/** @} */ 
