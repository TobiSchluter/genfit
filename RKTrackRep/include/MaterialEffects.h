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

#include "RKTools.h"
#include "GFPointPath.h"
#include <GFAbsMaterialInterface.h>

#include <iostream>
#include <vector>

#include <TObject.h>
#include <TVector3.h>


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

  //! set the material interface here. Material interface classes must be derived from GFAbsMaterialInterface.
  void init(GFAbsMaterialInterface* matIfc);

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


  //! Calculates energy loss in the traveled path, optional calculation of noise matrix
  double effects(const std::vector<GFPointPath>& points,
                 const double& mom,
                 const int& pdg,
                 double& xx0,
                 double* noise7x7 = NULL,
                 const double* jacobian7x7 = NULL,
                 const TVector3* directionBefore = NULL,
                 const TVector3* directionAfter = NULL);

  //! Returns maximum length so that a specified momentum loss will not be exceeded
  /**  The stepper returns the maximum length that the particle may travel, so that a specified relative momentum loss will not be exceeded,
   *   or the next material boundary is reached. The material crossed are stored together with their stepsizes.
  */
  double stepper(const RKTrackRep* rep,
                 M1x7& state7,
                 double sMax, // maximum step. unsigned!
                 const double& mom, // momentum
                 double& relMomLoss, // relative momloss for the step will be added
                 const int& pdg,
                 bool varField = true);


private:
  //! sets fcharge, fmass and calculates fbeta, fgamma, fgammasquare;
  void getParticleParameters(double mom);

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
  /**  This function first calcuates a MSC variance based on the current material and step length
   * 2 different formulas for the MSC variance are implemeted. One can select the formula via "setMscModel".
   * With the MSC variance and the current direction of the track a full 7D noise matrix is calculated.
   * This noise matrix is the additional noise at the end of fStep in the 7D globa cooridnate system
   * taking even the (co)variances of the position coordinates into account.
   * 
    */
  void noiseCoulomb(const double& mom,
                    double* noise,
                    const TVector3& direction) const;

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

  GFAbsMaterialInterface* fMaterialInterface;

public:
  ClassDef(GFMaterialEffects, 4);

};

#endif

/** @} */
