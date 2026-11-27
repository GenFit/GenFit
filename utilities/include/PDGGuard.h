#pragma once

#include <mutex>

namespace genfit {

/**
 * Shared guard for ROOT's TDatabasePDG singleton.
 *
 * TDatabasePDG::Instance() creates the database lazily and GetParticle() may add
 * entries lazily, so both must not run concurrently from several fit threads. All
 * lookups in genfit go through this one mutex (a per-file mutex would not protect
 * the ROOT singleton itself).
 */
inline std::mutex& pdgSingletonMutex() {
  static std::mutex m;
  return m;
}

} // namespace genfit
