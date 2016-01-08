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


#ifndef genfit_FieldManager_h
#define genfit_FieldManager_h

#include "AbsBField.h"
#include "IO.h"

#include <stdexcept>
#include <string>

#define CACHE

namespace genfit {

#ifdef CACHE
/**
 * @brief Cache B field at a position. Used by FieldManager.
 */
struct fieldCache {
  double posX; double posY; double posZ;
  double Bx; double By; double Bz;
};
#endif


/** @brief Singleton which provides access to magnetic field maps.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 */
class FieldManager {

 public:

  AbsBField* getField(){
    checkInitialized();
    return field_;
  }

  //! This does NOT use the cache!
  TVector3 getFieldVal(const TVector3& position){
    checkInitialized();
    return field_->get(position);
  }

#ifdef CACHE
  void getFieldVal(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz);
#else
  inline void getFieldVal(const double& posX, const double& posY, const double& posZ, double& Bx, double& By, double& Bz) {
    checkInitialized();
    return field_->get(posX, posY, posZ, Bx, By, Bz);
  }
#endif

  //! set the magnetic field here. Magnetic field classes must be derived from AbsBField.
  void init(AbsBField* b) {
    field_=b;
  }

  bool isInitialized() { return field_ != NULL; }

  void checkInitialized() {
    if(! isInitialized()){
      errorOut << "FieldManager hasn't been initialized with a correct AbsBField pointer!" << std::endl;
      std::string msg("FieldManager hasn't been initialized with a correct AbsBField pointer!");
      std::runtime_error err(msg);
      throw err;
    }
  }

  static void checkInstanciated() {
    if(instance_==NULL){
      errorOut << "FieldManager hasn't been instantiated yet, call getInstance() and init() before getFieldVal()!" << std::endl;
      std::string msg("FieldManager hasn't been instantiated yet, call getInstance() and init() before getFieldVal()!");
      std::runtime_error err(msg);
      throw err;
    }
  }

#ifdef CACHE
  //! Cache last lookup positions, and use stored field values if a lookup at (almost) the same position is done.
  void useCache(bool opt = true, unsigned int nBuckets = 8);
#else
  void useCache(bool opt = true, unsigned int nBuckets = 8) {
    std::cerr << "genfit::FieldManager::useCache() - FieldManager is compiled w/o CACHE, no caching will be done!" << std::endl;
  }
#endif

  //! Get singleton instance.
  static FieldManager* getInstance(){
    if(instance_ == NULL) {
      instance_ = new FieldManager();
    }
    return instance_;
  }


 private:

  FieldManager() {}
#ifdef CACHE
  ~FieldManager() { delete cache_; }
#else
  ~FieldManager() { }
#endif
  static FieldManager* instance_;
  static AbsBField* field_;

#ifdef CACHE
  static bool useCache_;
  static unsigned int n_buckets_;
  static fieldCache* cache_;
#endif

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_FieldManager_h
