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

#pragma once

#include "AbsBField.h"
#include "IO.h"

#include <mutex>
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
  ROOT::Math::XYZVector getFieldVal(const ROOT::Math::XYZVector& position){
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
    recordConfig();
  }

  void destruct() {
    {
      std::lock_guard<std::mutex> lock(configMutex_);
      delete configTemplate_;
      configTemplate_ = nullptr;
    }
    delete instance_;
    instance_ = nullptr;
  }

  bool isInitialized() { return field_ != nullptr; }

  void checkInitialized() {
    if(! isInitialized()){
      errorOut << "FieldManager hasn't been initialized with a correct AbsBField pointer!" << std::endl;
      std::string msg("FieldManager hasn't been initialized with a correct AbsBField pointer!");
      std::runtime_error err(msg);
      throw err;
    }
  }

  static void checkInstanciated() {
    if(instance_==nullptr){
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

  //! Get the calling thread's instance. One instance per thread: the field pointer and
  //! the cache settings are configuration, shared through a template; the lookup cache
  //! itself is per instance (per thread) and therefore never shared.
  static FieldManager* getInstance(){
    if(instance_ == nullptr) {
      std::lock_guard<std::mutex> lock(configMutex_);
      if (configTemplate_ == nullptr) configTemplate_ = new FieldManager();
      instance_ = new FieldManager(*configTemplate_);
    }
    return instance_;
  }


 private:

  FieldManager() {}
  //! Copies the configuration and allocates a fresh lookup cache.
  FieldManager(const FieldManager& cfg);
#ifdef CACHE
  ~FieldManager() { delete cache_; }
#else
  ~FieldManager() { }
#endif

  //! Copy the configuration fields to the template (thread-safe).
  void recordConfig();

  // One instance per thread: the field-lookup cache is per instance (per thread) and
  // therefore never shared between threads. The configuration set through init() and
  // useCache() is kept in a template which every instance created later copies, so the
  // public API (getInstance()) is unchanged.
  static thread_local FieldManager* instance_;
  static FieldManager* configTemplate_;
  static std::mutex configMutex_;
  AbsBField* field_ = nullptr;

#ifdef CACHE
  bool useCache_ = false;
  unsigned int n_buckets_ = 8;
  fieldCache* cache_ = nullptr;
  // cache sweep indices (previously function-local statics, i.e. shared state)
  int last_read_i_ = 0;
  int last_written_i_ = 0;
#endif

};

} /* End of namespace genfit */
/** @} */
