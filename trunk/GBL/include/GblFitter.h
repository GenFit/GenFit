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

#ifndef GblFitter_H
#define GblFitter_H

#include "GblTrajectory.h"
#include "AbsFitter.h"
#include "AbsTrackRep.h"
#include "GblFitterInfo.h"
#include "GblFitStatus.h"
#include "GblTrackSegmentController.h"

#include <map>
#include <iostream>

#include <TMatrixD.h>
#include <assert.h>
#include <sstream>

#include <TMath.h>
#include <TVector3.h>


namespace genfit {
  
  class GblTrackSegmentController;
  
  /** @brief Generic GBL implementation
   * 
   * The interface class to GBL track fit
   *
   */
  class GblFitter : public AbsFitter {
    
  private:
    GblFitter(const GblFitter&);
    GblFitter& operator=(GblFitter const&);
    
    std::string m_gblInternalIterations;
    bool m_enableScatterers;
    bool m_enableIntermediateScatterer;
    unsigned int m_externalIterations;
    unsigned int m_recalcJacobians;
    
    // Minimum scattering sigma (will be squared and inverted...)
    double scatEpsilon;
    GblTrackSegmentController* m_segmentController;
    
  public:
    
    /**
     * Default (and only) constructor
     */
    GblFitter() : AbsFitter(), m_gblInternalIterations(""), m_enableScatterers(true), m_enableIntermediateScatterer(true), m_externalIterations(1), m_recalcJacobians(0), scatEpsilon(1.e-8), m_segmentController(NULL) {;}
    
    /**
     * Destructor
     */
    virtual ~GblFitter();    

    /**
     * @brief Set options of the fitter/GBL
     * 
     * @return void
     * @param internalIterations GBL down-weighting in iterations. One letter (T,H,C) per iteration.
     *                           Seems "HH" is resonable for outliers. Default "" is OK.
     *                           Separate by comma for each external iteration
     *                           (usually last), e.g., for 3 iterations: ",,HH" to down-weight at last one
     * @param enableScatterers If false, no scatterers will be added to GBL trajectory
     * @param enableIntermediateScatterer True to simulate thick sctatterers by two thin scatterers
     *                                    1st at detector plane and intermediate between each two planes
     * @param externalIterations Sets number of times GblTrajectory.fit(...) will be called in processTrackWithRep(...).
     *                           Called external iterations. State is updated each time with GBL results.
     *                           If =0, GblFitterInfos will filled by reference states and GblFitStatus in track able
     *                           to construct simple GBL trajectory (for fit by GBL of output to Mille binary) / get list of points
     *                                    1st at detector plane and intermediate between each two planes
     * @param recalcJacobians Number of iteration up to which Jacobians should be recalculated / planes/meas updated after the fit.
     *                        0 = do not recalculate Jacobians. 1 = recalculate after first GBL fit. 2 = after 1st and 2nd GBL fit etc.
     */
    void setOptions(std::string internalIterations = "", bool enableScatterers = true, bool enableIntermediateScatterer = true, unsigned int externalIterations = 1, unsigned int recalcJacobians = 1) {
      m_externalIterations = externalIterations;
      m_gblInternalIterations = internalIterations;
      m_recalcJacobians = recalcJacobians;
      if (!enableScatterers)
        enableIntermediateScatterer = false;
      m_enableScatterers = enableScatterers;
      m_enableIntermediateScatterer = enableIntermediateScatterer;
    }

    
    /**
     * @brief Set multiple scattering options of the fitter/GBL
     * 
     * @return void
     * @param enableScatterers If false, no scatterers will be added to GBL trajectory
     * @param enableIntermediateScatterer True to simulate thick sctatterers by two thin scatterers
     *                                    1st at detector plane and intermediate between each two planes
     */
    void setMSOptions(bool enableScatterers = true, bool enableIntermediateScatterer = true) {
      if (!enableScatterers)
        enableIntermediateScatterer = false;
      m_enableScatterers = enableScatterers;
      m_enableIntermediateScatterer = enableIntermediateScatterer;
    }
    
    /**
     * @brief Evaluates moments of radiation length distribution from list of
     * material steps and computes parameters describing a corresponding thick scatterer.
     *
     * Based on input from Claus Kleinwort (DESY),
     * adapted for continuous material distribution represented by
     * a sum of step functions. Returned thick scatterer can be represented by two GBL scattering points.
     * Calculates variance of theta from total sum of radiation lengths
     * instead of summimg squares of individual deflection angle variances.
     *
     * @param length returned: Length of the track
     * @param theta returned: Variation of distribution of deflection angle
     * @param s returned: First moment of material scattering distribution
     * @param ds returned: Second moment (variance) of material scattering distribution
     * @param p Particle momentum magnitude (GeV/c)
     * @param mass Mass of particle (GeV/c/c)
     * @param steps Vector of material steps from (RKTrackRep) extrapolation
     * @return void
     */
    void getScattererFromMatList(double& length,
                                 double& theta, double& s, double& ds,
                                 const double p, const double mass, const double charge,
                                 const std::vector<genfit::MatStep>& steps) const;
    
    /**
     * Performs fit on a Track.
     * Hit resorting currently supported (use only if necessary /wire chamber/ ... will 
     * extrapolate along whole track to sort the hits).
     */
    void processTrackWithRep(Track* trk, const AbsTrackRep* rep, bool resortHits = false);
    
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
     * @return const std::vector< gbl::GblPoint, std::allocator >&
     */
    std::vector<gbl::GblPoint> collectGblPoints(genfit::Track* trk, const genfit::AbsTrackRep* rep);
    
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
    
    void setTrackSegmentController(GblTrackSegmentController* controler);
    
    
  public:
    
    ClassDef(GblFitter, 2)
    
  };
  
}  /* End of namespace genfit */
/** @} */

#endif // GblFitter_H

