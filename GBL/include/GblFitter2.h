/* Copyright 2013
 *   Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *   This file is part of GENFIT.
 *
 *   GENFIT is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   GENFIT is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */
/** @addtogroup genfit
 * @{
 */

#ifndef GblFitter2_H
#define GblFitter2_H

#include "GblTrajectory.h"
#include "AbsFitter.h"
#include "AbsTrackRep.h"
#include "GblFitterInfo2.h"
#include "GblFitStatus.h"

#include <map>
#include <iostream>

#include <TMatrixD.h>
#include <assert.h>
#include <sstream>

#include <TMath.h>
#include <TVector3.h>


namespace genfit {
    
  /** @brief Generic GBL implementation
   * 
   * The interface class to GBLv3 track fit
   *
   */
  class GblFitter2 : public AbsFitter {
    
  private:
    GblFitter2(const GblFitter2&);
    GblFitter2& operator=(GblFitter2 const&);
    
    std::string m_gblInternalIterations;
    unsigned int m_externalIterations;
    unsigned int m_resolveAmbiguties;
   
  public:
    
    /**
     * Default (and only) constructor
     */
    GblFitter2() : AbsFitter(), m_gblInternalIterations(""), m_externalIterations(1), m_resolveAmbiguties(0) {;}
    
    /**
     * Destructor
     */
    virtual ~GblFitter2();    

    /**
     * @brief Set options of the fitter/GBL
     * 
     * @return void
     * @param internalIterations GBL down-weighting in iterations. One letter (T,H,C) per iteration.
     *                           Seems "HH" is resonable for outliers. Default "" is OK.
     *                           Separate by comma for each external iteration
     *                           (usually last), e.g., for 3 iterations: ",,HH" to down-weight at last one
     * @param externalIterations Sets number of times GblTrajectory.fit(...) will be called in processTrackWithRep(...).
     *                           Called external iterations. State is updated each time with GBL results.
     *                           If =0, GblFitterInfos will filled by reference states and GblFitStatus in track able
     *                           to construct simple GBL trajectory (for fit by GBL of output to Mille binary) / get list of points
     *                                    1st at detector plane and intermediate between each two planes
     * @param recalcJacobians Number of iteration up to which Jacobians should be recalculated / planes/meas updated after the fit.
     *                        0 = do not recalculate Jacobians. 1 = recalculate after first GBL fit. 2 = after 1st and 2nd GBL fit etc.
     */
    void setOptions(const std::string &internalIterations = "", unsigned int externalIterations = 1, unsigned int resolveAmbiguties = 1) {
      m_externalIterations = externalIterations;
      m_gblInternalIterations = internalIterations;
      m_resolveAmbiguties = resolveAmbiguties;
    }
    
    /**
     * Performs fit on a Track.
     * Hit resorting currently supported (use only if necessary /wire chamber/ ... will 
     * extrapolate along whole track to sort the hits).
     */
    void processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits = false) override;
    
    /**
     * @brief Propagate seed, populate track with scatterers
     * and GblFitterInfos with reference state set
     * 
     * @param trk Track to attach with infos at given rep
     * @param rep TrackRep to which fitter info shall be attached
     * @return Length of track from extrapolations
     */
    double constructGblInfo(Track* trk, const AbsTrackRep* rep);
    
    /**
     * @brief Update down weights from trajectory fit. 
     *
     * @param traj The fitted GblTrajectory
     * @param trk The track with fitter infos from whose points traj was created 
     * @param rep The representation to which this fit status belong
     */
    void updateGblDownweights(gbl::GblTrajectory& traj, genfit::Track* trk, const genfit::AbsTrackRep* rep);        
    
    /**
     * @brief Populate all fitter infos in track for rep with
     * results of trajectory fit. 
     *
     * Updates also seed state in track (from forward prediction
     * at first point)
     * 
     * (The trajectory can only be cut before track end,
     * cannot have missing points in between (if valid))
     * 
     * 
     * 
     * TODO ??
     * Re-construct all points in GblFitterInfos (updated) and collect them
     * in fit status
     * 
     * 
     * @param traj The fitted GblTrajectory
     * @param trk The track with fitter infos from whose points traj was created 
     * @param rep The representation to which this fit status belong
     */
    void updateGblInfo(gbl::GblTrajectory& traj, genfit::Track* trk, const genfit::AbsTrackRep* rep);
        
    /**
     * @brief Constructs all GBL points and returns them in vector
     * for trajectory construction
     * 
     * @param trk The track
     * @param rep Representation for which to clean data
     * @param allowAmbiguities Allow ambiguities (to be resolved by GBL)
     * @return const std::vector< gbl::GblPoint, std::allocator >&
     */
    std::vector<gbl::GblPoint> collectGblPoints(genfit::Track* trk, const genfit::AbsTrackRep* rep, bool allowAmbiguities);
    
    /**
     * @brief Remove all previous gbl fitter data from track
     * Also removes trackpoints without measurement
     * 
     * @param trk The track
     * @param rep Representation for which to clean data
     * @return void
     */
    void cleanGblInfo(Track* trk, const AbsTrackRep* rep) const;
    
    /**
     * @brief Sort hits in track by arc-len using extrapolation
     * 
     * @param trk The track
     * @param rep The track representation
     * @return void
     */
    void sortHits(Track* trk, const AbsTrackRep* rep) const;
        
    
  public:
    
    ClassDef(GblFitter2, 1)
    
  };
  
}  /* End of namespace genfit */
/** @} */

#endif // GblFitter2_H

