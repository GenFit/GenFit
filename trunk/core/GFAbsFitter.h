/* Copyright 2012, Ludwig-Maximilians Universität München

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

#ifndef GFABSFITTER_H
#define GFABSFITTER_H

#include <algorithm>

class GFTrack;

/** @brief Common interface of track fitters.
 *
 */
class GFAbsFitter
{
 public:
  // Constructors/Destructors ---------
  GFAbsFitter() : fBlowUpFactor(500.) {}
  virtual ~GFAbsFitter() {}

  // Operators
  /** @brief Operator for use with STL.
   *
   * This operator allows to use the std::foreach algorithm with an
   * STL container of GFTrack* objects.
   */
  void operator()(GFTrack* track) { this->processTrack(track); }
  
  /** @brief Operator for use with STL.
   *
   * This operator allows to use the std::foreach algorithm with an
   * STL container of GFTrack* objects.
   */
  void operator()(std::pair<int,GFTrack*> tr) { this->processTrack(tr.second); }

  // Operations ----------------------

  /** @brief Performs fit on a GFTrack.
   *
   * The hits are processed in the order in which they are stored in the GFTrack
   * object. Sorting of hits in space has to be done before!
   */
  virtual void processTrack(GFTrack* trk) = 0;

  /** @brief Set the blowup factor (see blowUpCovs() )
   */
  void setBlowUpFactor(double f){ fBlowUpFactor=f; }
  
  /** @brief Get the blowup factor (see blowUpCovs() )
   */
  double getBlowUpFactor() const { return fBlowUpFactor; }

  /** @brief this is needed to blow up the covariance matrix before a fitting pass
   * drops off-diagonal elements and blows up diagonal by blowUpFactor
   */
  void blowUpCovs(GFTrack* trk) const;

 private:
  double fBlowUpFactor;
};

#endif

/** @} */
