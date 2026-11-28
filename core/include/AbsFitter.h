/* Copyright 2013, Ludwig-Maximilians Universität München,
   Authors: Tobias Schlüter & Johannes Rauch

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

#include <mutex>

namespace genfit {

class Track;
class AbsTrackRep;

/**
 * @brief Abstract base class for fitters.
 */
class AbsFitter {
 public:
  AbsFitter() : debugLvl_(0) {}
  virtual ~AbsFitter() {}

  //! Per-object operation lock: serializes concurrent calls on this fitter.
  std::recursive_mutex& usageMutex() const { return usageMutex_; }

  /**
   * Process Track with one AbsTrackRep of the Track. Optionally resort the hits if necessary (and supported by the fitter)
   */
  virtual void processTrackWithRep(Track*, const AbsTrackRep*, bool resortHits = false) = 0;

  /**
   * Process all reps. Start with the cardinalRep and resort the hits if necessary (and supported by the fitter)
   */
  void processTrack(Track*, bool resortHits = false);

  virtual void setDebugLvl(unsigned int lvl = 1) {debugLvl_ = lvl;}


 protected:

  unsigned int debugLvl_;

 private:

  mutable std::recursive_mutex usageMutex_;

};

}  /* End of namespace genfit */
/** @} */
