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
#include <iostream>
#include <assert.h>

#include <TObject.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TDatabasePDG.h>

#include <cmath>

/** @brief Track candidate -- a list of cluster indices
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
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
 * However this information is not autmatically used in genfit!!!
 * But a pointer to a GFTrackCand can be passed to the RKTrackRep to make use of this information.
 *
 * @sa RecoHitFactory
 */

class GFTrackCand : public TObject {
public:

  // Constructors/Destructors ---------
  GFTrackCand();
  ~GFTrackCand();

  //  /** @brief Initializing constructor
  //   *
  //   * @param curv Curvature from prefit. There is no stringent definition what
  //   * this parameter means at the moment.
  //   * @param dip Dip angle from prefit. There is no stringent definition what
  //   * this parameter means at the moment.
  //   * @param inv Dummy paramter. Has been used to mark inverted tracks
  //   * in the past.
  //   * @param detIDs collection of detector IDs. Each detector ID needs
  //   * a corresponding GFRecoHitProducer. See RecoHitFactory for details.
  //   * @param hitIDs collection of hit indices.
  //   */
  //  GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs);
  //  /* @brief same as previous ctor, but with ordering parameters */
  //  GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs, std::vector<double> rhos);

  /* @brief == operator does not check for rho */
  friend bool operator== (const GFTrackCand& lhs, const GFTrackCand& rhs);

  // Accessors -----------------------
  /** @brief Get detector ID and cluster index (hitId) for hit number i
   */
  void getHit(unsigned int i,
              unsigned int& detId,
              unsigned int& hitId) const {
    assert(i < getNHits());
    detId = fDetId.at(i); hitId = fHitId.at(i);
  }
  /** @brief Get detector ID and cluster index (hitId) for
   * hit number i with ordering parameter rho
   */
  void getHit(unsigned int i,
              unsigned int& detId,
              unsigned int& hitId,
              double& rho) const {
    assert(i < getNHits());
    detId = fDetId.at(i); hitId = fHitId.at(i);
    rho = fRho.at(i);
  }
  /** @brief Get detector ID and cluster index (hitId) for
   * hit number i with plane id
   */
  void getHitWithPlane(unsigned int i,
                       unsigned int& detId,
                       unsigned int& hitId,
                       unsigned int& planeId) const {
    assert(i < getNHits());
    detId = fDetId.at(i); hitId = fHitId.at(i);
    planeId = fPlaneId.at(i);
  }

  unsigned int getNHits() const {return fDetId.size();}
  double getCurv() const {return fCurv;}
  double getDip() const {return fDip;}
  //bool inverted() const {return fInv;} //nobody seems to use it so I commented it out
  std::vector<unsigned int> getHitIDs(int detId = -1) const;
  std::vector<unsigned int> GetHitIDs(int detId = -1) const;
  std::vector<unsigned int> getDetIDs() const {return fDetId;}
  std::vector<unsigned int> GetDetIDs() const { std::cerr << "the method GFTrackCand::GetDetIDs() is deprecated. Use GFTrackCand::getDetIDs() instead\n"; return fDetId;}
  std::vector<double>       getRhos() const {return fRho;}
  std::vector<double>       GetRhos() const {
    std::cerr << "the method GFTrackCand::GetRhos() is deprecated. Use GFTrackCand::getRhos() instead\n";
    return fRho;
  }
  std::set<unsigned int> getUniqueDetIDs() const {
    std::set<unsigned int> retVal;
    for (unsigned int i = 0; i < fDetId.size(); ++i) {
      retVal.insert(fDetId.at(i));
    }
    return retVal;
  }
  /** @brief get the MCT track id, for MC simulations - def. value -1
   */
  int getMcTrackId() const {return fMcTrackId;}
  /** @brief get the seed value for track: pos */
  TVector3 getPosSeed() const {
    std::cerr << "the method GFTrackCand::getPosSeed() is deprecated. Use GFTrackCand::getStateSeed() instead\n";
    TVector3 posSeed(fState6D[0][0], fState6D[1][0], fState6D[2][0]);
    return posSeed;
  }
  /** @brief get the seed value for track: direction */
  TVector3 getDirSeed() const {
    std::cerr << "the method GFTrackCand::getDirSeed() is deprecated. Use GFTrackCand::getStateSeed() instead\n";
    TVector3 dirSeed(fState6D[3][0], fState6D[4][0], fState6D[5][0]);
    return dirSeed;
  }
  /** @brief get the seed value for track: qoverp */
  double getQoverPseed() const {
    std::cerr << "the method GFTrackCand::getQoverPseed() is deprecated. Use GFTrackCand::getStateSeed() and/or getChargeSeed() instead\n";
    double p = std::sqrt(fState6D[3][0] * fState6D[3][0] + fState6D[4][0] * fState6D[4][0] + fState6D[5][0] * fState6D[5][0]);
    return fQ / p;
  }
  //  /** @brief get the seed value for track: error on pos (standard deviation)*/
  //  TVector3 getPosError() const {return fPosError;}
  //  /** @brief get the seed value for track: error on direction (standard deviation)*/
  //  TVector3 getDirError() const {return fDirError;}
  /** returns the 6D seed state; should be in global coordinates */
  TMatrixD getStateSeed() const {
    return fState6D;
  }
  /** returns the 6D covariance matrix of the seed state; should be in global coordinates */
  TMatrixD getCovSeed() const {
    return fCov6D;
  }
  double getChargeSeed()const {
    return fQ;
  }
  /** @brief get the PDG code*/
  int getPdgCode() const {return fPdg;}
  // Modifiers -----------------------
  void addHit(unsigned int detId, unsigned int hitId, double rho = 0., unsigned int planeId = 0);
  void setCurv(double c) {
    std::cerr << "the method GFTrackCand::setCurv is deprecated. Use GFTrackCand::set6DSeed instead\n";
    fCurv = c;
  }
  void setDip(double d) {
    std::cerr << "the method GFTrackCand::setDip is deprecated. Use GFTrackCand::set6DSeed instead\n";
    fDip = d;
  }
  //void setInverted(bool f=true) {fInv=f;} //nobody seems to use it so I commented it out
  /** @brief set the MCT track id, for MC simulations
   */
  void setMcTrackId(int i) {fMcTrackId = i;}
  /** @brief Test if hit already is part of this track candidate
   */
  bool hitInTrack(unsigned int detId, unsigned int hitId) const;
  /** @brief set the seed values for track: pos, direction, q/p
   */
  //  void setTrackSeed(const TVector3& pos,const TVector3& direction,const double qop){
  //    fPosSeed=pos;fDirSeed=direction;fQoverpSeed=qop;
  //  }
  /** @brief set the seed values for track: pos, momentum, pdgCode, pos error, momentum error (errors are optional and will be set to 1,1,1 if not given)
   */
  void setComplTrackSeed(const TVector3& pos, const TVector3& mom, const int pdgCode, TVector3 posError = TVector3(1.0, 1.0, 1.0), TVector3 dirError = TVector3(1.0, 1.0, 1.0));
  /** @brief set a particle hypothesis in form of a PDG code this will also set the charge attribute
   */
  void setPdgCode(int pdgCode) {
    fPdg = pdgCode;
    TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(fPdg);
    fQ = part->Charge() / (3.);
  }
  void append(const GFTrackCand&);

  /** @brief sort the hits that were already added to the trackCand using the rho parameter.
   * After this function was called rho will determine the order of propagation not the order of the addHit calls
   */
  void sortHits();

  void sortHits(std::vector<unsigned int> indices);

  // Operations ----------------------
  void reset();
  void Print(const Option_t* = "") const ;
  /** @brief sets the state to seed the track fitting. State has to be a TMatrixD(6,1). First 3 elements are the staring postion second 3 elements the starting momentum. Everything in global coordinates
   * charge is the charge hypotheses of the particle charge
   * ATTENTION: If you set the cov6D covariance matrix of the state remember that there ar VARIANCES not STANDARD DEVIATIONS on the diagonal
   */
  void set6DSeed(const TMatrixD& state6D, const double charge, TMatrixD cov6D = -1.0 * TMatrixD(TMatrixD::kUnit, TMatrixD(6, 6))) {
    fQ = charge;
    fState6D = state6D;
    fCov6D = cov6D;
  }
  /** @brief This function works the same as set6DSeed but instead of a charge hypothesis you can set a pdg code which will set the charge automatically
   * ATTENTION: If you set the cov6D covariance matrix of the state remember that there are VARIANCES not standard deviations on the diagonal
   */
  void set6DSeedAndPdgCode(const TMatrixD& state6D, const int pdgCode, TMatrixD cov6D = -1.0 * TMatrixD(TMatrixD::kUnit, TMatrixD(6, 6))) {
    setPdgCode(pdgCode);
    fState6D = state6D;
    fCov6D = cov6D;
  }
//  void setCurvDipPosSeed(const double curvature, const double dipAngle, const TVector3& pos){
//
//  };
//  /** set the errors  */
//  void setCurvDipPosError(const double curvatureError, const double dipAngleError, const TVector3& posError){};

private:

  // Private Data Members ------------
  std::vector<unsigned int> fDetId;
  std::vector<unsigned int> fHitId;
  std::vector<unsigned int> fPlaneId;
  std::vector<double>       fRho;

  double fCurv; // curvature from pattern reco
  double fDip;  // dip angle from pattern reco
  //bool fInv;  // true if inverted track //nobody seems to use it so I commented it out

//  TVector3 fPosSeed;  //seed value for the track: pos
//  TVector3 fDirSeed;  //direction
//  TVector3 fPosError;  //error on position seed given as a standard deviation
//  TVector3 fDirError;  //error on direction seed given as a standard deviation
//  double fQoverpSeed; //q/p

  int fMcTrackId; //if MC simulation, store the mct track id here
  int fPdg; // particle data groupe's id for a particle

  TMatrixD fState6D; /** global 6D position plus momentum state */
  TMatrixD fCov6D; /** global 6D position plus momentum covariance matrix */
  double fQ; /** the charge of the particle in units of elementary charge */


  // Private Methods -----------------

public:
  ClassDef(GFTrackCand, 8)
};

#endif

/** @} */
