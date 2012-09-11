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
#include "assert.h"

#include "TObject.h"
#include "TVector3.h"

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

  /** @brief Initializing constructor
   *
   * @param curv Curvature from prefit. There is no stringent definition what
   * this parameter means at the moment.
   * @param dip Dip angle from prefit. There is no stringent definition what
   * this parameter means at the moment.
   * @param inv Dummy paramter. Has been used to mark inverted tracks 
   * in the past.
   * @param detIDs collection of detector IDs. Each detector ID needs
   * a corresponding GFRecoHitProducer. See RecoHitFactory for details.
   * @param hitIDs collection of hit indices. 
   */
  GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs);
  /* @brief same as previous ctor, but with ordering parameters */
  GFTrackCand(double curv, double dip, double inv, std::vector<unsigned int> detIDs, std::vector<unsigned int> hitIDs, std::vector<double> rhos);

  /* @brief == operator does not check for rho */
  friend bool operator== (const GFTrackCand& lhs, const GFTrackCand& rhs);

  // Accessors -----------------------
  /** @brief Get detector ID and cluster index (hitId) for hit number i 
   */
  void getHit(unsigned int i, 
	      unsigned int& detId,
	      unsigned int& hitId) const {
	assert(i<getNHits());
	detId=fDetId.at(i);hitId=fHitId.at(i);
  }
  /** @brief Get detector ID and cluster index (hitId) for 
   * hit number i with ordering parameter rho 
   */
  void getHit(unsigned int i, 
	      unsigned int& detId,
	      unsigned int& hitId,
	      double &rho) const {
	assert(i<getNHits());
	detId=fDetId.at(i);hitId=fHitId.at(i);
	rho=fRho.at(i);
  }
  /** @brief Get detector ID and cluster index (hitId) for 
   * hit number i with plane id
   */
  void getHitWithPlane(unsigned int i, 
	      unsigned int& detId,
	      unsigned int& hitId,
	      unsigned int& planeId) const {
	assert(i<getNHits());
	detId=fDetId.at(i);hitId=fHitId.at(i);
	planeId=fPlaneId.at(i);
  }

  unsigned int getNHits() const {return fDetId.size();}
  double getCurv() const {return fCurv;}
  double getDip() const {return fDip;}
  bool inverted() const {return fInv;}
  std::vector<unsigned int> GetHitIDs(int detId=-1) const;
  std::vector<unsigned int> GetDetIDs() const {return fDetId;}
  std::vector<double>       GetRhos() const {return fRho;}
  std::set<unsigned int> GetUniqueDetIDs() const {
    std::set<unsigned int> retVal;
    for(unsigned int i=0;i<fDetId.size();++i){
      retVal.insert(fDetId.at(i));
    }
    return retVal;
  }
  /** @brief get the MCT track id, for MC simulations - def. value -1
   */
  int getMcTrackId() const {return fMcTrackId;}
  /** @brief get the seed value for track: pos */
  TVector3 getPosSeed() const {return fPosSeed;}
  /** @brief get the seed value for track: direction */
  TVector3 getDirSeed() const {return fDirSeed;}
  /** @brief get the seed value for track: qoverp */
  double getQoverPseed() const {return fQoverpSeed;}
  /** @brief get the seed value for track: error on pos (standard deviation)*/
  TVector3 getPosError() const {return fPosError;}
  /** @brief get the seed value for track: error on direction (standard deviation)*/
  TVector3 getDirError() const {return fDirError;}
  /** @brief get the PDG code*/
  int getPdgCode() const {return fPdg;}
  // Modifiers -----------------------
  void addHit(unsigned int detId, unsigned int hitId, double rho=0., unsigned int planeId=0);
  void setCurv(double c){fCurv=c;}
  void setDip(double d){fDip=d;}
  void setInverted(bool f=true) {fInv=f;}
  /** @brief set the MCT track id, for MC simulations
   */
  void setMcTrackId(int i){fMcTrackId=i;}
  /** @brief Test if hit already is part of this track candidate
   */
  bool HitInTrack(unsigned int detId, unsigned int hitId) const;
  /** @brief set the seed values for track: pos, direction, q/p
   */
  void setTrackSeed(const TVector3& pos,const TVector3& direction,const double qop){
    fPosSeed=pos;fDirSeed=direction;fQoverpSeed=qop;
  }
  /** @brief set the seed values for track: pos, momentum, pdgCode, pos error, momentum error (errors are optional and will be set to 1,1,1 if not given)
   */
  void setComplTrackSeed(const TVector3& pos,const TVector3& mom, const int pdgCode, TVector3 posError = TVector3(1.0,1.0,1.0), TVector3 dirError = TVector3(1.0,1.0,1.0));
  /** @brief set a particle hypothesis in form of a PDG code
   */
  void setPdgCode(int pdgCode){fPdg=pdgCode;}
  void append(const GFTrackCand&);

  /** @brief sort the hits that were already added to the trackCand using the rho parameter.
   * After this function was called rho will determine the order of propagation not the order of the addHit calls
   */
  void sortHits();

  void sortHits(std::vector<unsigned int> indices);

  // Operations ----------------------
  void reset();
  void Print(const Option_t* = "") const ;

private:

  // Private Data Members ------------
  std::vector<unsigned int> fDetId;
  std::vector<unsigned int> fHitId;
  std::vector<unsigned int> fPlaneId;
  std::vector<double>       fRho;

  double fCurv; // curvature from pattern reco
  double fDip;  // dip angle from pattern reco
  bool fInv;  // true if inverted track

  TVector3 fPosSeed;  //seed value for the track: pos
  TVector3 fDirSeed;  //direction
  TVector3 fPosError;  //error on position seed given as a standard deviation
  TVector3 fDirError;  //error on direction seed given as a standard deviation
  double fQoverpSeed; //q/p
  int fMcTrackId; //if MC simulation, store the mct track id here
  int fPdg; // particle data groupe's id for a particle


  // Private Methods -----------------

public:
  ClassDef(GFTrackCand,7)
};

#endif

/** @} */
