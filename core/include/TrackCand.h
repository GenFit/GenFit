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
 * @{ */

#ifndef genfit_TrackCand_h
#define genfit_TrackCand_h

#include "TrackCandHit.h"

#include <vector>
#include <set>
#include <assert.h>

#include <TObject.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDatabasePDG.h>

#include <cmath>


namespace genfit {

/** @brief Track candidate -- seed values and indices
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Moritz Nadler (maintainer during 2012)
 *
 * The main task of the TrackCand object is to store a list of indices to
 * cluster objects. Each cluster in the Track is identified by it's
 * detector ID and it's index in the corresponding TClonesArray.
 * Also there is a ordering parameter to order hits.
 * Optionally, plane indices for the hits can be stored (most importantly
 * for fitting with the Daf).
 * This information is used by the RecoHitFactory to automatically load
 * RecoHits into a Track. Through this it is possible to define Tracks over
 * an arbitrary number of different detectors.
 *
 * In addition TrackCand offers members to store starting values for the fit.
 * The starting values (seeds) for the fit are stored as a 6D state (x,y,z,px,py,pz) and its
 * corresponding 6x6 covariance matrix. All seed getter and setter manipulate these two members
 * but the user can chose using TVector3 or TMatrixD to get/set the seed state.
 * However this information is not automatically used in genfit.
 * But a pointer to a TrackCand can be passed to the a RKTrackRep constructor
 * to make use of this information without manually extracting it from the TrackCand object.
 *
 * @sa RecoHitFactory
 */
class TrackCand : public TObject {


 public:


  // Constructors/Destructors ---------
  TrackCand();
  ~TrackCand();

  //! copy constructor
  TrackCand( const TrackCand& other );
  //! assignment operator
  TrackCand& operator=(TrackCand other);
  void swap(TrackCand& other); // nothrow

  //! == operator checks equality of TrackCandHits. Does not check for sorting parameters.
  friend bool operator== (const TrackCand& lhs, const TrackCand& rhs);
  friend bool operator!= (const TrackCand& lhs, const TrackCand& rhs) {return !(lhs == rhs);}

  static bool compareTrackCandHits(const TrackCandHit* lhs, const TrackCandHit* rhs) {return (*lhs < *rhs);} // operator< defined in TrackCandHit.h

  // Accessors -----------------------
  TrackCandHit* getHit(int i) const;

  //!Get detector Id and hit Id for hit number i
  void getHit(int i, int& detId, int& hitId) const;

  //! Get detector Id, hit Id and sorting parameter for hit number i
  void getHit(int i, int& detId, int& hitId, double& sortingParameter) const;

  //! Get detector Id, hit Id and plane id for hit number i
  void getHitWithPlane(int i, int& detId, int& hitId, int& planeId) const;

  unsigned int getNHits() const {return hits_.size();}

  /**
   * @brief Get hit ids of from a specific detector.
   *
   * DetId -1 gives hitIds of hits with default detId -1. The default argument -2 gives hit Ids of all hits.
   */
  std::vector<int>    getHitIDs(int detId = -2) const;

  //! Get detector IDs of all hits
  std::vector<int>    getDetIDs() const;
  //! Get sorting parameterts of all hits
  std::vector<double> getSortingParameters() const;
  std::set<int>       getUniqueDetIDs() const;

  //! Get the MCT track id, for MC simulations - default value = -1
  int getMcTrackId() const {return mcTrackId_;}

  //! Get the time at which the seed state is defined
  double getTimeSeed() const { return time_; }

  /** @brief get the seed value for track: pos. Identical to the first 3 components of getStateSeed*/
  TVector3 getPosSeed() const {return TVector3(state6D_(0), state6D_(1), state6D_(2));}

  /** @brief get the seed value for track: mom. Identical to the last 3 components of getStateSeed*/
  TVector3 getMomSeed() const {return TVector3(state6D_(3), state6D_(4), state6D_(5));}

  /** @brief get the covariance matrix seed (6D).  */
  const TMatrixDSym& getCovSeed() const {return cov6D_;}

  //! Returns the 6D seed state; should be in global coordinates.
  const TVectorD& getStateSeed() const {return state6D_;}

  double getChargeSeed() const {return q_;}

  //! Get the PDG code
  int getPdgCode() const {return pdg_;}

  //! Is there a hit with detId and hitId in the TrackCand?
  bool hitInTrack(int detId, int hitId) const;

  // Modifiers -----------------------

  void addHit(int detId, int hitId, int planeId = -1, double sortingParameter = 0);

  void addHit(TrackCandHit* hit) {hits_.push_back(hit);}

  //! Set the MCT track id, for MC simulations
  void setMcTrackId(int i) {mcTrackId_ = i;}

  //! Set a particle hypothesis in form of a PDG code. This will also set the charge attribute
  void setPdgCode(int pdgCode);

  //! Clone the TrackCandHit objects from the other TrackCand and append them to this TrackCand
  void append(const TrackCand&);

  //! Sort the hits that were already added to the trackCand using the sorting parameters.
  void sortHits();

  void sortHits(const std::vector<unsigned int>& indices);

  // Operations ----------------------
  //! Delete and clear the TrackCandHits
  void reset();

  //! Write the content of all private attributes to the terminal
  void Print(const Option_t* = "") const ;

  //! Set the time at which the seed is defined
  void setTimeSeed(double time) { time_ = time; }

  /** @brief set the covariance matrix seed (6D).  */
  void setCovSeed(const TMatrixDSym& cov6D) {cov6D_ = cov6D; /* always 6D, no need to resize */}

  /** @brief sets the state to seed the track fitting. State has to be a TVectorD(6). First 3 elements are the staring postion second 3 elements the starting momentum. Everything in global coordinates
   * charge is the charge hypotheses of the particle charge
   */
  void set6DSeed(const TVectorD& state6D, const double charge);

  /** @brief This function works the same as set6DSeed but instead of a charge hypothesis you can set a pdg code which will set the charge automatically
   */
  void set6DSeedAndPdgCode(const TVectorD& state6D, const int pdgCode);

  /** @brief sets the state to seed the track fitting. State has to be a TVector3 for position and a TVector3 for momentum. Everything in global coordinates
   * charge is the charge hypotheses of the particle charge
   */
  void setPosMomSeed(const TVector3& pos, const TVector3& mom, const double charge);

  /** @brief This function works the same as setPosMomSeed but instead of a charge hypothesis you can set a pdg code which will set the charge automatically
   */
  void setPosMomSeedAndPdgCode(const TVector3& pos, const TVector3& mom, const int pdgCode);

  /** @brief sets the state to seed the track fitting and its
     time. State has to be a TVectorD(6). First 3 elements are the
     staring postion second 3 elements the starting
     momentum. Everything in global coordinates charge is the charge
     hypotheses of the particle charge.
   */
  void setTime6DSeed(double time, const TVectorD& state6D, const double charge);

  /** @brief This function works the same as set6DSeed but instead of
      a charge hypothesis you can set a pdg code which will set the
      charge automatically.
   */
  void setTime6DSeedAndPdgCode(double time, const TVectorD& state6D, const int pdgCode);

  /** @brief sets the state to seed the track fitting and its time. State has to be
     a TVector3 for position and a TVector3 for momentum. Everything
     in global coordinates charge is the charge hypotheses of the
     particle charge.
   */
  void setTimePosMomSeed(double time, const TVector3& pos, const TVector3& mom,
			 const double charge);

  /** @brief This function works the same as setPosMomSeed but instead
      of a charge hypothesis you can set a pdg code which will set the
      charge automatically.
   */
  void setTimePosMomSeedAndPdgCode(double time, const TVector3& pos,
				   const TVector3& mom, const int pdgCode);


 private:

  // Private Data Members ------------
  std::vector<TrackCandHit*> hits_; //->

  int mcTrackId_; /**< if MC simulation, store the mc track id here */
  int pdg_; /**< particle data groupe's id for a particle*/

  double time_; /**< Time at which the seed is given */
  TVectorD state6D_; /**< global 6D position plus momentum state */
  TMatrixDSym cov6D_; /**< global 6D position plus momentum state */
  double q_; /**< the charge of the particle in units of elementary charge */


 public:

  ClassDef(TrackCand,2)
  // Version history:
  // ver 2: keep track of time in state (schema evolution rule added).
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackCand_h
