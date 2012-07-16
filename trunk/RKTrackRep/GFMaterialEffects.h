/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef GFMATERIALEFFECTS_H
#define GFMATERIALEFFECTS_H

#include<iostream>
#include "RKTools.h"
#include "TObject.h"
#include <vector>
#include "TVector3.h"
#include "TMatrixT.h"
#include "GFPointPath.h"

/** @brief  Contains stepper and energy loss/noise matrix calculation
 *
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

class GFMaterialEffects : public TObject {
private:
  GFMaterialEffects();
  virtual ~GFMaterialEffects();

  static GFMaterialEffects* finstance;

public:
  static GFMaterialEffects* getInstance();
  static void destruct();

  void setNoEffects(bool opt = true) {fNoEffects = opt;}

  void setEnergyLossBetheBloch(bool opt = true) {fEnergyLossBetheBloch = opt; fNoEffects = false;}
  void setNoiseBetheBloch(bool opt = true) {fNoiseBetheBloch = opt; fNoEffects = false;}
  void setNoiseCoulomb(bool opt = true) {fNoiseCoulomb = opt; fNoEffects = false;}
  void setEnergyLossBrems(bool opt = true) {fEnergyLossBrems = opt; fNoEffects = false;}
  void setNoiseBrems(bool opt = true) {fNoiseBrems = opt; fNoEffects = false;}

  /** @brief Select the multiple scattering model that will be used during track fit.
   *  At the moment two model are available GEANE and Highland. GEANE is the model was was present in Genfit first.
   *  Note that using this function has no effect if setNoiseCoulomb(false) is set.
   */
  void setMscModel(const std::string& modelName);


  //! Calculates energy loss in the travelled path, optional calculation of noise matrix
  double effects(const std::vector<GFPointPath>& points,
                 const double& mom,
                 const int& pdg,
                 double& xx0,
                 const bool& doNoise = false,
                 double* noise7x7 = NULL,
                 const double* jacobian7x7 = NULL,
                 const TVector3* directionBefore = NULL,
                 const TVector3* directionAfter = NULL);

  //! Returns maximum length so that a specified momentum loss will not be exceeded
  /**  The stepper returns the maximum length that the particle may travel, so that a specified relative momentum loss will not be exceeded,
   *   or the next material boundary is reached. The material crossed are stored together with their stepsizes.
  */
  double stepper(const double& maxStep, // maximum step. unsigned!
                 const double& maxAngleStep, // maximum step due to curvature. unsigned!
                 const double& posx,
                 const double& posy,
                 const double& posz,
                 const double& dirx,
                 const double& diry,
                 const double& dirz,
                 const double& mom, // momentum
                 double& relMomLoss, // relative momloss for the step will be added
                 const int& pdg);


private:
  //! sets fcharge, fmass and calculates fbeta, fgamma, fgammasquare;
  void getParticleParameters(double mom);

  //! sets fmatDensity, fmatZ, fmatA, fradiationLength, fmEE;
  void getMaterialParameters(TGeoMaterial* mat);

  //! Returns energy loss
  /**  Uses Bethe Bloch formula to calculate energy loss.
    *  Calcuates and sets fdedx which needed also for noiseBetheBloch.
    *  Therefore it is not a const function!
    *
  */
  double energyLossBetheBloch(const double& mom);

  //! calculation of energy loss straggeling
  /**  For the energy loss straggeling, different formulas are used for different regions:
    *  - Vavilov-Gaussian regime
    *  - Urban/Landau approximation
    *  - truncated Landau distribution
    *  - Urban model
    *
    *  Needs fdedx, which is calculated in energyLossBetheBloch, so it has to be called afterwards!
    */
  void noiseBetheBloch(const double& mom,
                       double* noise) const;

  //! calculation of multiple scattering
  /**  With the calculated multiple scattering angle, two noise matrices are calculated:
    *  - with respect to #directionBefore: #noiseBefore, which is then propagated with the jacobian
    *  - with respect to #directionAfter: #noiseAfter
    *  The mean value of these two matrices is added to the noise matrix #noise.
    *  This method gives better results than either calculating only noiseBefore or noiseAfter.
    *  \n\n
    *  This is a detailed description of the mathematics involved:
    *  \n
    *
    *  We define a local coordinate system cs' with initial momentum in z-direction:
    *  \f[
    *  \left(\begin{array}{c}
    *  x'\\
    *  y'\\
    *  z'\\
    *  a_{x}'\\
    *  a_{y}'\\
    *  a_{z}'\\
    *  \frac{q}{p}'\end{array}\right)=\left(\begin{array}{ccccccc}
    *  \cos\psi & \sin\psi & 0\\
    *  -\cos\vartheta\sin\psi & \cos\vartheta\cos\psi & \sin\psi\\
    *  \sin\vartheta\sin\psi & -\sin\vartheta\cos\psi & \cos\vartheta\\
    *   &  &  & \cos\psi & \sin\psi & 0\\
    *   &  &  & -\cos\vartheta\sin\psi & \cos\vartheta\cos\psi & \sin\psi\\
    *   &  &  & \sin\vartheta\sin\psi & -\sin\vartheta\cos\psi & \cos\vartheta\\
    *   &  &  &  &  &  & 1\end{array}\right)\left(\begin{array}{c}
    *  x\\
    *  y\\
    *  z\\
    *  a_{x}\\
    *  a_{y}\\
    *  a_{z}\\
    *  \frac{q}{p}\end{array}\right)=R^{T}\left(\begin{array}{c}
    *  x\\
    *  y\\
    *  z\\
    *  a_{x}\\
    *  a_{y}\\
    *  a_{z}\\
    *  \frac{q}{p}\end{array}\right)\f]
    *
    *  \n
    *
    * Now the global coordinate system cs is:\n
    *  \f[
    *  \left(\begin{array}{c}
    *  x\\
    *  y\\
    *  z\\
    *  a_{x}\\
    *  a_{y}\\
    *  a_{z}\\
    *  \frac{q}{p}\end{array}\right)=\left(\begin{array}{ccccccc}
    *  \cos\psi & -\cos\vartheta\sin\psi & \sin\vartheta\sin\psi\\
    *  \sin\psi & \cos\vartheta\cos\psi & -\sin\vartheta\cos\psi\\
    *  0 & \sin\psi & \cos\vartheta\\
    *   &  &  & \cos\psi & -\cos\vartheta\sin\psi & \sin\vartheta\sin\psi\\
    *   &  &  & \sin\psi & \cos\vartheta\cos\psi & -\sin\vartheta\cos\psi\\
    *   &  &  & 0 & \sin\psi & \cos\vartheta\\
    *   &  &  &  &  &  & 1\end{array}\right)\left(\begin{array}{c}
    *  x'\\
    *  y'\\
    *  z'\\
    *  a_{x}'\\
    *  a_{y}'\\
    *  a_{z}'\\
    *  \frac{q}{p}'\end{array}\right)=R\left(\begin{array}{c}
    *  x'\\
    *  y'\\
    *  z'\\
    *  a_{x}'\\
    *  a_{y}'\\
    *  a_{z}'\\
    *  \frac{q}{p}'\end{array}\right) \f]
    *
    *  \n
    *
    *  with the Euler angles
    *
    * \f{eqnarray*}
    * \psi & = &
    * \begin{cases}
    *   \begin{cases}
    *     \frac{\pi}{2} & a_{x} \geq 0 \\
    *     \frac{3\pi}{2} & a_{x} < 0
    *   \end{cases}
    *   & a_{y}=0 \mbox{ resp. } |a_{y}|<10^{-14} \\
    *   - \arctan \frac{a_{x}}{a_{y}} & a_{y} < 0 \\
    *   \pi - \arctan \frac{a_{x}}{a_{y}} & a_{y} > 0
    * \end{cases} \\
    * \vartheta & = & \arccos a_{z}
    * \f}
    *
    * \n
    * \f$M\f$ is the multiple scattering error-matrix in the \f$\theta\f$ coordinate system.
    * \f$\theta_{1/2}=0\f$ are the multiple scattering angles. There is only an error in direction
    * (which later leads to an error in position, when \f$N_{before}\f$ is propagated with \f$T\f$).
    * \f[
    * M=\left(\begin{array}{cc}
    * \sigma^{2} & 0\\
    * 0 & \sigma^{2}\end{array}\right)\f]
    *
    * \n
    * This translates into the noise matrix \f$\overline{M}\f$ in the local cs' via jacobian J.
    * The mean scattering angle is always 0, this means that \f$\theta_{1/2}=0\f$):
    *
    * \f{eqnarray*}
    * x' & = & x'\\
    * y' & = & y'\\
    * z' & = & z'\\
    * a_{x}' & = & \sin\theta_{1}\\
    * a_{y}' & = & \sin\theta_{2}\\
    * a_{z}' & = & \sqrt{1-\sin^{2}\theta_{1}-\sin^{2}\theta_{2}}\\
    * \frac{q}{p}' & = & \frac{q}{p}'\f}
    *
    *
    * \f[
    * M=\left(\begin{array}{cc}
    * \sigma^{2} & 0\\
    * 0 & \sigma^{2}\end{array}\right)\f]
    *
    * \n
    *
    * \f[
    * J=\frac{\partial\left(x',y',z',a_{x}',a_{y}',a_{z}',\frac{q}{p}'\right)}{\partial\left(\theta_{1},\theta_{2}\right)}\f]
    *
    * \n
    *
    * \f[
    * J=\left(\begin{array}{cc}
    * 0 & 0\\
    * 0 & 0\\
    * 0 & 0\\
    * \cos\theta_{1} & 0\\
    * 0 & \cos\theta_{2}\\
    * -\frac{\cos\theta_{1}\sin\theta_{1}}{\sqrt{\cos^{2}\theta_{1}-\sin^{2}\theta_{2}}} & -\frac{\cos\theta_{2}\sin\theta_{2}}{\sqrt{\cos^{2}\theta_{1}-\sin^{2}\theta_{2}}}\\
    * 0 & 0\end{array}\right) \overset{\theta_{1/2}=0}{=} \left(\begin{array}{cc}
    * 0 & 0\\
    * 0 & 0\\
    * 0 & 0\\
    * 1 & 0\\
    * 0 & 1\\
    * 0 & 0\\
    * 0 & 0\end{array}\right)\f]
    * \n
    *
    * \f[
    * \overline{M}=J\: M\: J^{T}=\left(\begin{array}{ccccccc}
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & \sigma^{2} & 0 & 0 & 0\\
    * 0 & 0 & 0 & 0 & \sigma^{2} & 0 & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\end{array}\right)\f]
    *
    *
    * The following transformation brings the noise matrix into the global coordinate system cs, resulting in  \f$N\f$ :
    *
    *\f[
    * N=R\overline{M}R^{T}=\sigma^{2}\left(\begin{array}{ccccccc}
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\\
    * 0 & 0 & 0 & \cos^{2}\psi+\cos^{2}\theta-\cos^{2}\theta\cos^{2}\psi & \cos\psi\sin\psi\sin^{2}\theta & -\cos\theta\sin\psi\sin\theta & 0\\
    * 0 & 0 & 0 & \cos\psi\sin\psi\sin^{2}\theta & \sin^{2}\psi+\cos^{2}\theta\cos^{2}\psi & \cos\theta\cos\psi\sin\theta & 0\\
    * 0 & 0 & 0 & -\cos\theta\sin\psi\sin\theta & \cos\theta\cos\psi\sin\theta & \sin^{2}\theta & 0\\
    * 0 & 0 & 0 & 0 & 0 & 0 & 0\end{array}\right)\f]
    *
    * \n
    *
    *
    * Now two \f$N\f$ are calculated: \f$N_{before}\f$ (at the start point, with initial direction #directionBefore)
    * and \f$N_{after}\f$ (at the final point, with direction #directionAfter).
    * \f$N_{before}\f$ is the propagated with the #jacobian of extrapolation \f$T\f$.
    * The new covariance matrix with multiple scattering in global cs is:
    *
    * \f{eqnarray*}
    * C_{new} & = C_{old} + 0.5 \cdot T N_{before} T^{T} + 0.5 \cdot N_{after} \f}
    *
    *
    * \n
    * \n
    * See also: Track following in dense media and inhomogeneous magnetic fields,
    * A.Fontana, P.Genova, L.Lavezzi and A.Rotondi;
    * PANDA Report PV/01-07
    * \n
    * \n
    */
  void noiseCoulomb(const double& mom,
                    double* noise,
                    const double* jacobian,
                    const TVector3* directionBefore,
                    const TVector3* directionAfter) const;

  //! Returns energy loss
  /** Can be called with any pdg, but only calculates energy loss for electrons and positrons (otherwise returns 0).
    * Uses a gaussian approximation (Bethe-Heitler formula with Migdal corrections).
    * For positrons the energy loss is weighed with a correction factor.
  */
  double energyLossBrems(const double& mom) const;

  //! calculation of energy loss straggeling
  /** Can be called with any pdg, but only calculates straggeling for electrons and positrons.
   *
   */
  void noiseBrems(const double& mom,
                  double* noise) const;


  bool fNoEffects;

  bool fEnergyLossBetheBloch;
  bool fNoiseBetheBloch;
  bool fNoiseCoulomb;
  bool fEnergyLossBrems;
  bool fNoiseBrems;

  const double me; // electron mass (GeV)

  double fstep; // stepsize

  // cached values for energy loss and noise calculations
  double fbeta;
  double fdedx;
  double fgamma;
  double fgammaSquare;

  double fmatDensity;
  double fmatZ;
  double fmatA;
  double fradiationLength;
  double fmEE; // mean excitation energy

  int fpdg;
  double fcharge;
  double fmass;

  int fMscModelCode; /// depending on this number a specific msc model is chosen in the noiseCoulomb function.

public:
  ClassDef(GFMaterialEffects, 2);

};

#endif

/** @} */
