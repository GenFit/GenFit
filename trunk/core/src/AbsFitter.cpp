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

#include "AbsFitter.h"
#include "Track.h"

namespace genfit {

void AbsFitter::processTrack(Track* tr, bool resortHits) {
  AbsTrackRep* cardRep = tr->getCardinalRep();
  // process cardinal rep first
  processTrackWithRep(tr, cardRep, resortHits);

  // now process rest of reps, but don't change sorting anymore!
  for (unsigned int i=0; i<tr->getNumReps(); ++i) {
    if (tr->getTrackRep(i) != cardRep)
      processTrackWithRep(tr, tr->getTrackRep(i), false);
  }

  // self check
  assert(tr->checkConsistency());
}

} /* End of namespace genfit */
