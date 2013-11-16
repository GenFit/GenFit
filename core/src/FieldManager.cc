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
#include "FieldManager.h"

#include <iostream>
#include <math.h>

namespace genfit {

FieldManager* FieldManager::instance_ = NULL;
AbsBField* FieldManager::field_ = NULL;
bool FieldManager::useCache_ = false;
unsigned int FieldManager::n_buckets_ = 8;
fieldCache* FieldManager::cache_ = NULL;

//#define DEBUG

void FieldManager::getFieldVal(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz){
  checkInitialized();

  if (useCache_) {

    // cache code copied from http://en.wikibooks.org/wiki/Optimizing_C%2B%2B/General_optimization_techniques/Memoization
    static int last_read_i = 0;
    static int last_written_i = 0;
    int i = last_read_i;

    static const double epsilon = 0.001;

    #ifdef DEBUG
    static int used = 0;
    static int notUsed = 0;
    #endif

    do {
      if (fabs(cache_[i].posX - posX) < epsilon &&
	  fabs(cache_[i].posY - posY) < epsilon &&
	  fabs(cache_[i].posZ - posZ) < epsilon) {
	Bx = cache_[i].Bx;
	By = cache_[i].By;
	Bz = cache_[i].Bz;
        #ifdef DEBUG
	++used;
	std::cout<<"used the cache! " << double(used)/(used + notUsed) << "\n";
        #endif
	return;
      }
      i = (i + 1) % n_buckets_;
    } while (i != last_read_i);

    last_read_i = last_written_i = (last_written_i + 1) % n_buckets_;

    cache_[last_written_i].posX = posX;
    cache_[last_written_i].posY = posY;
    cache_[last_written_i].posZ = posZ;

    field_->get(posX, posY, posZ, cache_[last_written_i].Bx, cache_[last_written_i].By, cache_[last_written_i].Bz);

    Bx = cache_[last_written_i].Bx;
    By = cache_[last_written_i].By;
    Bz = cache_[last_written_i].Bz;
    #ifdef DEBUG
    ++notUsed;
    std::cout<<"did NOT use the cache! \n";
    #endif
    return;

  }
  else
    return field_->get(posX, posY, posZ, Bx, By, Bz);

}


void FieldManager::useCache(bool opt, unsigned int nBuckets) {
  useCache_ = opt;
  n_buckets_ = nBuckets;

  if (useCache_) {
    cache_ = new fieldCache[n_buckets_];
    for (size_t i = 0; i < n_buckets_; ++i) {
      // Should be safe to initialize with values in Andromeda
      cache_[i].posX = cache_[i].posY = cache_[i].posZ = 2.4e24 / sqrt(3);
      cache_[i].Bx = cache_[i].By = cache_[i].Bz = 1e30;      
    }      
  }
}


} /* End of namespace genfit */
