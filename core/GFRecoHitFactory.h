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
 * @{ */

/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 */

#ifndef GFRECOHITFACTORY_H
#define GFRECOHITFACTORY_H

#include<vector>
#include<map>

#include "GFRecoHitProducer.h"
#include "GFTrackCand.h"

class GFAbsRecoHit;


/** @brief Factory object to create RecoHits from digitized and clustered data
 *
 * The GFRecoHitFactory is used to automatically fill Track objects with 
 * hit data. For each detector type that is used, one GFRecoHitProducer 
 * has to be registered in the factory. The factory can the use the index
 * information from a GFTrackCand object to load the indexed hits into
 * the Track. 
 *
 * @sa GFAbsRecoHitProducer
 * @sa GFTrackCand
 */
class GFRecoHitFactory{
private:
  std::map<int,GFAbsRecoHitProducer*> fHitProdMap;


public:
  GFRecoHitFactory();
  virtual ~GFRecoHitFactory();

  /** @brief Register a producer module to the factory
   *
   * For each type of hit a separate producer is needed. The type of hit
   * is identified by the detector ID (detID). This index corresponds to the
   * detector ID that is stored in the GFTrackCand object
   */
  void addProducer(int detID, GFAbsRecoHitProducer* hitProd);

  /** @brief Clear all hit producers
   */
  void clear();

  /** @brief Create a RecoHit
   *
   * RecoHits have to implement a Constructor which takes the cluster object
   * from which the RecoHit is build as the only parameter. 
   * See GFAbsRecoHitProducer for details
   */
  GFAbsRecoHit*              createOne (int detID,int index);

  /** @brief Creat a collection of RecoHits
   *
   * This is the standard way to prepare the hit collection for a Track. The
   * resulting collection can contain hits from several detectors. The order
   * of the hits is the same as in the GFTrackCand. It is assumed that this order
   * is already along the track.
   *
   * RecoHits have to implement a constructor which takes the cluster object
   * from which the RecoHit is build as the only parameter. 
   * See GFAbsRecoHitProducer for details
   */
  std::vector<GFAbsRecoHit*> createMany(const GFTrackCand& cand);
					 

};


#endif

/** @} */
