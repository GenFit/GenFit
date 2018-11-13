/* Copyright 2008-2014, Technische Universitaet Muenchen,
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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_MaterialEffects_h
#define genfit_MaterialEffects_h

#include "RKTools.h"
#include "AbsMaterialInterface.h"

#include <iostream>
#include <vector>

#include <TVector3.h>


namespace genfit {

/** @brief Stepper and energy loss/noise matrix calculation
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, author)
 *
 *  It provides functionality to limit the stepsize of an extrapolation in order not to
 *  exceed a specified maximum momentum loss. After propagation, the energy loss
 *  for the given length and (optionally) the noise matrix can be calculated.
 *  You have to set which energy-loss and noise mechanisms you want to use.
 *  At the moment, per default all energy loss and noise options are ON.
 */
class MaterialEffects {

 private:

  MaterialEffects();
  virtual ~MaterialEffects();

  static MaterialEffects* instance_;


public:

  static MaterialEffects* getInstance();
  static void destruct();

  //! set the material interface here. Material interface classes must be derived from AbsMaterialInterface.
  void init(AbsMaterialInterface* matIfc);
  bool isInitialized() { return materialInterface_ != nullptr; }

  void setNoEffects(bool opt = true) {noEffects_ = opt;}

  void setEnergyLossBetheBloch(bool opt = true) {energyLossBetheBloch_ = opt; noEffects_ = false;}
  void setNoiseBetheBloch(bool opt = true) {noiseBetheBloch_ = opt; noEffects_ = false;}
  void setNoiseCoulomb(bool opt = true) {noiseCoulomb_ = opt; noEffects_ = false;}
  void setEnergyLossBrems(bool opt = true) {energyLossBrems_ = opt; noEffects_ = false;}
  void setNoiseBrems(bool opt = true) {noiseBrems_ = opt; noEffects_ = false;}
  void ignoreBoundariesBetweenEqualMaterials(bool opt = true) {ignoreBoundariesBetweenEqualMaterials_ = opt;}

  /** @brief Select the multiple scattering model that will be used during track fit.
   *
   *  At the moment two model are available GEANE and Highland. GEANE is the model was was present in Genfit first.
   *  Note that using this function has no effect if setNoiseCoulomb(false) is set.
   */
  void setMscModel(const std::string& modelName);


  //! Calculates energy loss in the traveled path, optional calculation of noise matrix
  double effects(const std::vector<RKStep>& steps,
                 int materialsFXStart,
                 int materialsFXStop,
                 const double& mom,
                 const int& pdg,
                 M7x7* noise = nullptr);

  /**  @brief Returns maximum length so that a specified momentum loss will not be exceeded.
   *
   * The stepper returns the maximum length that the particle may travel, so that a specified relative momentum loss will not be exceeded,
   * or the next material boundary is reached. The material crossed are stored together with their stepsizes.
  */
  void stepper(const RKTrackRep* rep,
               M1x7& state7,
               const double& mom, // momentum
               double& relMomLoss, // relative momloss for the step will be added
               const int& pdg,
               Material& currentMaterial,
               StepLimits& limits,
               bool varField = true);

  void setDebugLvl(unsigned int lvl = 1);


  void drawdEdx(int pdg = 11);

 private:

  //! sets charge_, mass_
  void getParticleParameters();

  void getMomGammaBeta(double Energy,
                       double& mom, double& gammaSquare, double& gamma, double& betaSquare) const;

  //! Returns momentum loss
  /**
   * Also sets dEdx_ and E_.
   */
  double momentumLoss(double stepSign, double mom, bool linear);

  //! Calculate dEdx for a given energy
  double dEdx(double Energy) const;


  //! Uses Bethe Bloch formula to calculate dEdx.
  double dEdxBetheBloch(double betaSquare, double gamma, double gammasquare) const;

  //! calculation of energy loss straggeling
  /**  For the energy loss straggeling, different formulas are used for different regions:
    *  - Vavilov-Gaussian regime
    *  - Urban/Landau approximation
    *  - truncated Landau distribution
    *  - Urban model
    *
    *  Needs dEdx_, which is calculated in momentumLoss, so it has to be called afterwards!
    */
  void noiseBetheBloch(M7x7& noise, double mom, double betaSquare, double gamma, double gammaSquare) const;

  //! calculation of multiple scattering
  /**  This function first calcuates a MSC variance based on the current material and step length
   * 2 different formulas for the MSC variance are implemeted. One can select the formula via "setMscModel".
   * With the MSC variance and the current direction of the track a full 7D noise matrix is calculated.
   * This noise matrix is the additional noise at the end of fStep in the 7D globa cooridnate system
   * taking even the (co)variances of the position coordinates into account.
   * 
    */
  void noiseCoulomb(M7x7& noise,
                    const M1x3& direction, double momSquare, double betaSquare) const;

  //! Returns dEdx
  /** Can be called with any pdg, but only calculates dEdx for electrons and positrons (otherwise returns 0).
    * Uses a gaussian approximation (Bethe-Heitler formula with Migdal corrections).
    * For positrons, dEdx is weighed with a correction factor.
  */
  double dEdxBrems(double mom) const;

  //! calculation of energy loss straggeling
  /** Can be called with any pdg, but only calculates straggeling for electrons and positrons.
   */
  void noiseBrems(M7x7& noise, double momSquare, double betaSquare) const;



  bool noEffects_;

  bool energyLossBetheBloch_;
  bool noiseBetheBloch_;
  bool noiseCoulomb_;
  bool energyLossBrems_;
  bool noiseBrems_;

  bool ignoreBoundariesBetweenEqualMaterials_;

  const double me_; // electron mass (GeV)

  double stepSize_; // stepsize

  // cached values for energy loss and noise calculations
  double dEdx_; // Runkge Kutta dEdx
  double E_; // Runge Kutta Energy
  double matDensity_;
  double matZ_;
  double matA_;
  double radiationLength_;
  double mEE_; // mean excitation energy

  int pdg_;
  int charge_;
  double mass_;

  int mscModelCode_; /// depending on this number a specific msc model is chosen in the noiseCoulomb function.

  AbsMaterialInterface* materialInterface_;

  unsigned int debugLvl_;

  // ClassDef(MaterialEffects, 1);

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_MaterialEffects_h
