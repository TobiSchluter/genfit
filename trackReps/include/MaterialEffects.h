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
#include "RKTrackRep.h"
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
  bool initTrack(double posX, double posY, double posZ,
                 double dirX, double dirY, double dirZ)
  {
    bool result = materialInterface_->initTrack(posX, posY, posZ, dirX, dirY, dirZ);
    return result;
  }
  double findNextBoundaryStraightLine(double sMax)
  {
    return materialInterface_->findNextBoundaryStraightLine(sMax);
  }

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
  double effects(const std::vector<RKStep>::const_iterator& materialsFXStart,
                 const std::vector<RKStep>::const_iterator& materialsFXStop,
                 const double& mom,
                 const int& pdg,
                 M7x7* noise = nullptr);

  /**  @brief Returns maximum length so that a specified momentum loss will not be exceeded.
   *
   * The stepper returns the maximum length that the particle may travel, so that a specified relative momentum loss will not be exceeded,
   * or the next material boundary is reached. The material crossed are stored together with their stepsizes.
  */
  void stepper(AbsTrackRep::internalExtrapolator& extrap,
               const double& mom, // momentum
               double& relMomLoss, // relative momloss for the step will be added
               const int& pdg,
               MaterialProperties& currentMaterial,
               StepLimits& limits);

  void setDebugLvl(unsigned int lvl = 1);


  void drawdEdx(int pdg = 11);

  //! Calculate dEdx for a given energy
  double dEdx(const MaterialProperties& material, double Energy) const;
  double d2EdxdE(const MaterialProperties& material, double Energy); // derivative

  void getMaterialProperties(MaterialProperties& material) { materialInterface_->getMaterialParameters(material); }

  //! sets charge_, mass_
  void getParticleParameters(int pdg);

private:


  void getMomGammaBeta(double Energy,
                       double& mom, double& gammaSquare, double& gamma, double& betaSquare) const;

  //! Returns momentum loss
  /**
   * Also sets dEdx_ and E_.
   */
  double momentumLoss(const MaterialProperties& material, double stepSign, double mom, bool linear);

  //! Uses Bethe Bloch formula to calculate dEdx.
  double dEdxBetheBloch(const MaterialProperties& material, double betaSquare, double gamma, double gammasquare) const;

  //! calculation of energy loss straggeling
  /**  For the energy loss straggeling, different formulas are used for different regions:
    *  - Vavilov-Gaussian regime
    *  - Urban/Landau approximation
    *  - truncated Landau distribution
    *  - Urban model
    *
    *  Needs dEdx_, which is calculated in momentumLoss, so it has to be called afterwards!
    */
  void noiseBetheBloch(M7x7& noise, const MaterialProperties& material, double mom, double betaSquare, double gamma, double gammaSquare) const;

  //! calculation of multiple scattering
  /**  This function first calcuates a MSC variance based on the current material and step length
   * 2 different formulas for the MSC variance are implemeted. One can select the formula via "setMscModel".
   * With the MSC variance and the current direction of the track a full 7D noise matrix is calculated.
   * This noise matrix is the additional noise at the end of fStep in the 7D globa cooridnate system
   * taking even the (co)variances of the position coordinates into account.
   * 
    */
  void noiseCoulomb(M7x7& noise, const MaterialProperties& material,
                    const M1x3& direction, double momSquare, double betaSquare) const;

  //! Returns dEdx
  /** Can be called with any pdg, but only calculates dEdx for electrons and positrons (otherwise returns 0).
    * Uses a gaussian approximation (Bethe-Heitler formula with Migdal corrections).
    * For positrons, dEdx is weighed with a correction factor.
  */
  double dEdxBrems(const MaterialProperties& material, double mom) const;

  //! calculation of energy loss straggeling
  /** Can be called with any pdg, but only calculates straggeling for electrons and positrons.
   */
  void noiseBrems(M7x7& noise, const MaterialProperties& material, double momSquare, double betaSquare) const;



  bool noEffects_;

  bool energyLossBetheBloch_;
  bool noiseBetheBloch_;
  bool noiseCoulomb_;
  bool energyLossBrems_;
  bool noiseBrems_;

  bool ignoreBoundariesBetweenEqualMaterials_;

  double stepSize_; // stepsize

  // cached values for energy loss and noise calculations
  double dEdx_; // Runkge Kutta dEdx
  double E_; // Runge Kutta Energy

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
