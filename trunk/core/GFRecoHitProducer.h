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

#ifndef GFRECOHITPRODUCER_H
#define GFRECOHITPRODUCER_H

#include<vector>
#include<map>
#include<assert.h>
#include<iostream>

#include "TClonesArray.h"

#include "GFException.h"

class GFAbsRecoHit;

/** @brief Abstract interface class for GFRecoHitProducer
 *
 * Defines the very basic interface of a producer.
 */
class GFAbsRecoHitProducer {
public:
  /** @brief Virtual abstract method to produce a RecoHit. 
   * Implemented in GFRecoHitProducer
   */
  virtual GFAbsRecoHit* produce(int index)=0;
  virtual ~GFAbsRecoHitProducer();
};


/** @brief Template class for a hit producer module
 *
 * A GFRecoHitProducer module is used by RecoHitFactory to create RecoHits for
 * one specific detector type. 
 *
 * It is assumed that each detector has as output of its digitization /
 * clustering some sort of cluster class which stores all information that
 * corresponds to a measured hit in that detector. The RecoHit producer 
 * converts this information into a class that can be handled by genfit.
 * This class is realized as a RecoHit (a class inherting from GFAbsRecoHit).
 *
 * In order to use the GFRecoHitProducer facility a
 * RecoHit has to implement a constructor which takes as an argument 
 * a pointer to the cluster class. This constructor serves as the initializing
 * constructor for the RecoHit.
 *
 * The GFRecoHitProducer will fetch the cluster objects from a TClonesArray and
 * use the initializing constructor to build the corresponding RecoHit. 
 *
 * @param hit_t template parameter specifying cluster class
 * @param recoHit_t template parameter specifying recoHit
 */
template <class hit_T,class recoHit_T>
class GFRecoHitProducer : public GFAbsRecoHitProducer {
 private:
  /** @brief pointer to array with cluster data */
  TClonesArray* hitArrayTClones;
  //std::vector<GFAbsRecoHit*>* hitArrayVector;
 public:

  /** @brief Constructor takes pointer to the cluster array */
  GFRecoHitProducer(TClonesArray*);
  //GFRecoHitProducer(std::vector<GFAbsRecoHit*>*);
  virtual ~GFRecoHitProducer();

  /** @brief Create a RecoHit from the cluster at position index 
   * in TClonesArray
   */
  virtual GFAbsRecoHit* produce(int index);	
};
/** @} */

template <class hit_T,class recoHit_T>
GFRecoHitProducer<hit_T,recoHit_T>::GFRecoHitProducer(TClonesArray* theArr) {
  hitArrayTClones = theArr;
  //hitArrayVector = NULL;
}
/*
template <class hit_T,class recoHit_T>
  GFRecoHitProducer<hit_T,recoHit_T>::GFRecoHitProducer(std::vector<GFAbsRecoHit*>* theArr) {
  hitArrayTClones = NULL;
  hitArrayVector = theArr;
}
*/

template <class hit_T,class recoHit_T>
GFRecoHitProducer<hit_T,recoHit_T>::~GFRecoHitProducer() {
  //we dont assume ownership over the hit arrays
}


template <class hit_T,class recoHit_T>
GFAbsRecoHit* GFRecoHitProducer<hit_T,recoHit_T>::produce(int index) {
  assert(hitArrayTClones!=NULL);
  //assert(hitArrayTClones!=NULL || hitArrayVector!=NULL);//at least one exists
  //assert(!(hitArrayTClones!=NULL && hitArrayVector!=NULL));//but not both
  //if(hitArrayTClones!=NULL){
    //the ROOT guys really use 0 and not NULL grrr...
  if(hitArrayTClones->At(index) == 0) {
    GFException e("In GFRecoHitProducer: index for hit in TClonesArray out of bounds",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  return ( new recoHit_T( (hit_T*) hitArrayTClones->At(index) ) );
  //}
  //else{//after assertions this is save: the hitArrayVector is good
  //  if(index >= hitArrayVector->size()) {
  //    GFException e("In GFRecoHitProducer: index for hit in std::vector out of bounds",__LINE__,__FILE__);
  //    e.setFatal();
  //    throw e;
  //  }
  //  return ( new recoHit_T( (hit_T*) hitArrayVector->at(index) ) );
  //}
}


#endif 


