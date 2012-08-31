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

#ifndef GFRECOHITIFC_H
#define GFRECOHITIFC_H


#include "TMatrixT.h"

#include "GFAbsRecoHit.h"
#include "GFDetPlane.h"


/** @brief RecoHit interface template class. Provides comfortable 
 * interface to create RecoHits
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * This class defines a comfortable interface to create hit classes in genfit.
 * It is a template class. The template parameter is used to specify a certain
 * basic type of hit:
 *  - GFRecoHitIfc<PlanarHitPolicy> a basic planar hit
 *  - GFRecoHitIfc<SpacepointHitPolicy> a basic space point hit
 *  - GFRecoHitIfc<WirepointHitPolicy> a basic hit on a wire
 *
 * To create a hit for a detector simply inherit from one of the options 
 * above and fill in your data. For details look at the respective 
 * HitPolicy documentations. You can also directly inherit from 
 * GFAbsRecoHit though this is not recommended. If a new hit geometry is needed
 * one should think about implementing a new HitPolicy for this type of hit.
 *
 * @sa PlanarHitPolicy
 * @sa SpacepointHitPolicy
 * @sa WirepointHitPolicy
 *
 * Implementation details: The actual implementations of the methods 
 * declared here can be found in the HitPolicy objects.
 */
template<class HitPolicy>
class GFRecoHitIfc : public GFAbsRecoHit{
 protected:
  HitPolicy fPolicy;

 public:
  
  /** @brief Constructor specifying dimension of hit coordinate vector
   */
  GFRecoHitIfc(int dim) : GFAbsRecoHit(dim){;}
  virtual ~GFRecoHitIfc(){;}

  /** @brief Returns the detector plane object for this hit and a given track
   * representation.
   *
   * The actutal code for this method depends on the hit geometry and is 
   * implemented in the HitPolicy
   * @sa PlanarHitPolicy
   * @sa SpacepointHitPolicy
   * @sa WirepointHitPolicy
   */
  virtual const GFDetPlane& getDetPlane(GFAbsTrackRep* rep){return fPolicy.detPlane(this,rep);}

  /** @brief Get measurement m,V and HMatrix
   *
   * Implementation in the HitPolicy
   */
  virtual void getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TMatrixT<double>& statePred,const TMatrixT<double>& covPred,TMatrixT<double>& m, TMatrixT<double>& V){
    static_cast<void>(rep);
    static_cast<void>(statePred);
    static_cast<void>(covPred);
    TMatrixT<double> mTemp = fPolicy.hitCoord(this,pl);
    TMatrixT<double> VTemp = fPolicy.hitCov(this,pl);
    m.ResizeTo(mTemp);
    V.ResizeTo(VTemp);
    m = mTemp;
    V = VTemp;
}

  const std::string& getPolicyName(){return fPolicy.getName();}

 public:
  ClassDef(GFRecoHitIfc,1);

};

#endif


/** @} */
