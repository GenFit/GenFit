/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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

#ifndef GFTRACKCAND_H
#define GFTRACKCAND_H

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

#include "GFTrackCandHit.h"

/** @brief Track candidate -- seed values and indices
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Moritz Nadler (maintainer during 2012)
 * The main task of the GFTrackCand object is to store a list of indices to
 * cluster objects. Each cluster in the Track is identified by it's
 * detector ID and it's index in the corresponding TClonesArray.
 * Also there is a ordering parameter rho, to order hits.
 * Optionally, plane indices for the hits can be stored (most importantly
 * for fitting with the GFDaf).
 * This information is used by the RecoHitFactory to automatically load
 * RecoHits into a Track. Through this it is possible to define Tracks over
 * an arbitrary number of different detectors.
 *
 * In addition GFTrackCand offers members to store starting values for the fit.
 * The starting values (seeds) for the fit are stored as a 6D state (x,y,z,px,py,pz) and its
 * corresponding 6x6 covariance matrix. All seed getter and setter manipulate these two members
 * but the user can chose using TVector3 or TMatrixD to get/set the seed state.
 * However this information is not automatically used in genfit.
 * But a pointer to a GFTrackCand can be passed to the a RKTrackRep constructor
 * to make use of this information without manually extracting it from the GFTRackCand object.
 *
 * @sa RecoHitFactory
 */

class GFTrackCand : public TObject {
 public:

  // Constructors/Destructors ---------
  GFTrackCand();
  ~GFTrackCand();

  //! copy constructor
  GFTrackCand( const GFTrackCand& other );
  //! assignment operator
  GFTrackCand& operator=( const GFTrackCand& other );

  /* @brief == operator checks equality of TrackCandHits. Does not check for rho. */
  friend bool operator== (const GFTrackCand& lhs, const GFTrackCand& rhs);
  friend bool operator!= (const GFTrackCand& lhs, const GFTrackCand& rhs) {
    return !(lhs == rhs);
  }

  static bool compareTrackCandHits(const GFTrackCandHit* lhs, const GFTrackCandHit* rhs){
    return (*lhs < *rhs); // operator< defined in GFTrackCandHit.h
  }

  // Accessors -----------------------
  GFTrackCandHit* getHit(unsigned int i) const {
    assert(i < getNHits());
    return fHits[i];
  }

  /** @brief Get detector Id and hit Id for hit number i
   */
  void getHit(unsigned int i,
              int& detId,
              int& hitId) const {
    assert(i < getNHits());
    detId = fHits[i]->getDetId();
    hitId = fHits[i]->getHitId();
  }

  /** @brief Get detector Id, hit Id and ordering
   * parameter rho for hit number i
   */
  void getHit(unsigned int i,
              int& detId,
              int& hitId,
              double& rho) const {
    assert(i < getNHits());
    detId = fHits[i]->getDetId();
    hitId = fHits[i]->getHitId();
    rho = fHits[i]->getRho();
  }

  /** @brief Get detector Id, hit Id and plane id
   * for hit number i
   */
  void getHitWithPlane(unsigned int i,
                       int& detId,
                       int& hitId,
                       int& planeId) const {
    assert(i < getNHits());
    detId = fHits[i]->getDetId();
    hitId = fHits[i]->getHitId();
    planeId = fHits[i]->getPlaneId();
  }

  unsigned int getNHits() const {return fHits.size();}

  /** @brief get hit ids of from a specific detector. DetId -1 gives hitIds of
   * hits with default detId -1. The default argument -2 gives hit Ids of all hits.
   */
  std::vector<int>    getHitIDs(int detId = -2) const;
  std::vector<int>    getDetIDs() const;
  std::vector<double> getRhos() const;
  std::set<int>       getUniqueDetIDs() const;

  /** @brief get the MCT track id, for MC simulations - default value -1
   */
  int getMcTrackId() const {return fMcTrackId;}

  /** @brief get the seed value for track: pos. Identical to the first 3 components of getStateSeed*/
  TVector3 getPosSeed() const {
    TVector3 posSeed(fState6D(0), fState6D(1), fState6D(2));
    return posSeed;
  }

  /** @brief get the seed value for track: mom. Identical to the last 3 components of getStateSeed*/
  TVector3 getMomSeed() const {
    TVector3 momSeed(fState6D(3), fState6D(4), fState6D(5));
    return momSeed;
  }

  /** returns the 6D seed state; should be in global coordinates */
  TVectorD getStateSeed() const {
    return fState6D;
  }

  /** returns the 6D covariance matrix of the seed state; should be in global coordinates */
  TMatrixDSym getCovSeed() const {
    return fCov6D;
  }

  double getChargeSeed()const {
    return fQ;
  }

  /** @brief get the PDG code*/
  int getPdgCode() const {return fPdg;}

  bool hitInTrack(int detId, int hitId) const;

  // Modifiers -----------------------
  void addHit(int detId,
              int hitId,
              int planeId = -1,
              double rho = 0);

  void addHit(GFTrackCandHit* hit) {fHits.push_back(hit);}

  /** @brief set the MCT track id, for MC simulations
   */
  void setMcTrackId(int i) {fMcTrackId = i;}

  /** @brief Test if hit already is part of this track candidate
   */

  /** @brief set a particle hypothesis in form of a PDG code. This will also set the charge attribute
   */
  void setPdgCode(int pdgCode) {
    fPdg = pdgCode;
    TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(fPdg);
    fQ = part->Charge() / (3.);
  }

  void append(const GFTrackCand&);

  /** @brief sort the hits that were already added to the trackCand using the rho parameter.
   */
  void sortHits();

  void sortHits(const std::vector<unsigned int>& indices);

  // Operations ----------------------
  /** @brief delete and clear the GFTrackCandHits
   */
  void reset();

  /** @brief write the content of all private attributes to the terminal
   */
  void Print(const Option_t* = "") const ;

  /** @brief sets the state to seed the track fitting. State has to be a TVectorD(6). First 3 elements are the staring postion second 3 elements the starting momentum. Everything in global coordinates
   * charge is the charge hypotheses of the particle charge
   * ATTENTION: If you set the cov6D covariance matrix of the state remember that there are VARIANCES not STANDARD DEVIATIONS on the diagonal
   */
  void set6DSeed(const TVectorD& state6D, const double charge, const TMatrixDSym& cov6D) {
    fQ = charge;
    fState6D = state6D;
    fCov6D = cov6D;
  }

  void set6DSeed(const TVectorD& state6D, const double charge) {
    fQ = charge;
    fState6D = state6D;
    fCov6D =  -1.0 * TMatrixDSym(TMatrixDSym::kUnit, TMatrixDSym(6));
  }

  /** @brief This function works the same as set6DSeed but instead of a charge hypothesis you can set a pdg code which will set the charge automatically
   * ATTENTION: If you set the cov6D covariance matrix of the state remember that there are VARIANCES not standard deviations on the diagonal
   */
  void set6DSeedAndPdgCode(const TVectorD& state6D, const int pdgCode, const TMatrixDSym& cov6D) {
    setPdgCode(pdgCode);
    fState6D = state6D;
    fCov6D = cov6D;
  }

  void set6DSeedAndPdgCode(const TVectorD& state6D, const int pdgCode) {
    setPdgCode(pdgCode);
    fState6D = state6D;
    fCov6D =  -1.0 * TMatrixDSym(TMatrixDSym::kUnit, TMatrixDSym(6));
  }

  /** @brief sets the state to seed the track fitting. State has to be a TVector3 for position and a TVector3 for momentum. Everything in global coordinates
   * charge is the charge hypotheses of the particle charge
   * ATTENTION: If you set the cov6D covariance matrix of the state remember that there are VARIANCES not STANDARD DEVIATIONS on the diagonal
   */
  void setPosMomSeed(const TVector3& pos, const TVector3& mom, const double charge, const TMatrixDSym& cov6D) {
    fQ = charge;
    fState6D[0] = pos[0];  fState6D[1] = pos[1];  fState6D[2] = pos[2];
    fState6D[3] = mom[0];  fState6D[4] = mom[1];  fState6D[5] = mom[2];
    fCov6D = cov6D;
  }

  void setPosMomSeed(const TVector3& pos, const TVector3& mom, const double charge) {
    fQ = charge;
    fState6D[0] = pos[0];  fState6D[1] = pos[1];  fState6D[2] = pos[2];
    fState6D[3] = mom[0];  fState6D[4] = mom[1];  fState6D[5] = mom[2];
    fCov6D =  -1.0 * TMatrixDSym(TMatrixDSym::kUnit, TMatrixDSym(6));
  }

  /** @brief This function works the same as setPosMomSeed but instead of a charge hypothesis you can set a pdg code which will set the charge automatically
   * ATTENTION: If you set the cov6D covariance matrix of the state remember that there are VARIANCES not standard deviations on the diagonal
   */
  void setPosMomSeedAndPdgCode(const TVector3& pos, const TVector3& mom, const int pdgCode, const TMatrixDSym& cov6D) {
    setPdgCode(pdgCode);
    fState6D[0] = pos[0];  fState6D[1] = pos[1];  fState6D[2] = pos[2];
    fState6D[3] = mom[0];  fState6D[4] = mom[1];  fState6D[5] = mom[2];
    fCov6D = cov6D;
  }

  void setPosMomSeedAndPdgCode(const TVector3& pos, const TVector3& mom, const int pdgCode) {
    setPdgCode(pdgCode);
    fState6D[0] = pos[0];  fState6D[1] = pos[1];  fState6D[2] = pos[2];
    fState6D[3] = mom[0];  fState6D[4] = mom[1];  fState6D[5] = mom[2];
    fCov6D = -1.0 * TMatrixDSym(TMatrixDSym::kUnit, TMatrixDSym(6));
  }

private:

  // Private Data Members ------------
  std::vector<GFTrackCandHit*> fHits; //->

  int fMcTrackId; /**< if MC simulation, store the mc track id here */
  int fPdg; /**< particle data groupe's id for a particle*/

  TVectorD fState6D; /**< global 6D position plus momentum state */
  TMatrixDSym fCov6D; /**< global 6D position plus momentum covariance matrix */
  double fQ; /**< the charge of the particle in units of elementary charge */


  // Private Methods -----------------

public:
  ClassDef(GFTrackCand, 11)
};

#endif

/** @} */
