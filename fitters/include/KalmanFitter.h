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

#ifndef genfit_KalmanFitter_h
#define genfit_KalmanFitter_h

#include "AbsKalmanFitter.h"

#include <memory>


namespace genfit {

class KalmanFitterInfo;
class MeasuredStateOnPlane;
class TrackPoint;

/**
 * @brief Simple Kalman filter implementation.
 */
class KalmanFitter : public AbsKalmanFitter {

 private:

  // These private functions are needed, otherwise strange things happen, no idea why!
  KalmanFitter(const KalmanFitter&);
  KalmanFitter& operator=(KalmanFitter const&);

 public:

  KalmanFitter(unsigned int maxIterations = 4, double deltaPval = 1e-3, double blowUpFactor = 1e3, bool squareRootFormalism = false)
    : AbsKalmanFitter(maxIterations, deltaPval, blowUpFactor), currentState_(nullptr),
      squareRootFormalism_(squareRootFormalism)
  {}

  ~KalmanFitter() {}

  //! Hit resorting currently NOT supported.
  void processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits = false);

  //! process only a part of the track. Can also be used to process the track only in backward direction.
  //! Does not alter the FitStatus and does not do multiple iterations.
  void processTrackPartially(Track* tr, const AbsTrackRep* rep, int startId = 0, int endId = -1);

  void useSquareRootFormalism(bool squareRootFormalism = true) {squareRootFormalism_ = squareRootFormalism;}

 private:
  bool fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int startId, int endId, int& nFailedHits);
  void processTrackPoint(TrackPoint* tp,
      const AbsTrackRep* rep, double& chi2, double& ndf, int direction);

  std::unique_ptr<MeasuredStateOnPlane> currentState_;

  bool squareRootFormalism_;

 public:
  ClassDef(KalmanFitter,1)

};

} /* End of namespace genfit */
/** @} */

#endif //genfit_KalmanFitter_h
