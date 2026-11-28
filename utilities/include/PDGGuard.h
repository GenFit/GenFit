#pragma once

#include <mutex>
#include <unordered_map>

#include <TDatabasePDG.h>

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

/**
 * Look up a particle in TDatabasePDG, caching the result per thread.
 *
 * The transport needs the same particle parameters on every Runge Kutta step. To
 * keep the shared mutex of pdgSingletonMutex() off that hot path, each thread
 * caches its lookup results here and the mutex is taken only on a cache miss, i.e.
 * at most once per thread and PDG code. Only successful lookups are cached: a
 * failed lookup is retried under the mutex so particles added to TDatabasePDG
 * later are still found.
 */
inline const TParticlePDG* pdgParticle(int pdgCode) {
  thread_local std::unordered_map<int, const TParticlePDG*> cache;
  const auto it = cache.find(pdgCode);
  if (it != cache.end()) return it->second;

  const TParticlePDG* particle;
  {
    std::lock_guard<std::mutex> lock(pdgSingletonMutex());
    particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
  }
  if (particle != nullptr) cache.emplace(pdgCode, particle);
  return particle;
}

} // namespace genfit
