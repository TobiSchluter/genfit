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

#include "MaterialEffects.h"
#include "Exception.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <TDatabasePDG.h>
#include <TMath.h>


//#define DEBUG

namespace genfit {

MaterialEffects* MaterialEffects::instance_ = nullptr;


MaterialEffects::MaterialEffects():
  noEffects_(false),
  energyLossBetheBloch_(true), noiseBetheBloch_(true),
  noiseCoulomb_(true),
  energyLossBrems_(true), noiseBrems_(true),
  me_(0.510998910E-3),
  stepSize_(0),
  beta_(0),
  dEdx_(0),
  gamma_(0),
  gammaSquare_(0),
  matDensity_(0),
  matZ_(0),
  matA_(0),
  radiationLength_(0),
  mEE_(0),
  pdg_(0),
  charge_(0),
  mass_(0),
  mscModelCode_(0),
  materialInterface_(nullptr)
{
}

MaterialEffects::~MaterialEffects()
{
  if (materialInterface_ != nullptr) delete materialInterface_;
}

MaterialEffects* MaterialEffects::getInstance()
{
  if (instance_ == nullptr) instance_ = new MaterialEffects();
  return instance_;
}

void MaterialEffects::destruct()
{
  if (instance_ != nullptr) {
    delete instance_;
    instance_ = nullptr;
  }
}

void MaterialEffects::init(AbsMaterialInterface* matIfc)
{
  if (materialInterface_ != nullptr) {
    std::string msg("MaterialEffects::initMaterialInterface(): Already initialized! ");
    std::runtime_error err(msg);
  }
  materialInterface_ = matIfc;
}



void MaterialEffects::setMscModel(const std::string& modelName)
{
  if (modelName == std::string("GEANE")) {
    mscModelCode_ = 0;
  } else if (modelName == std::string("Highland")) {
    mscModelCode_ = 1;
  } else {// throw exception
    std::string errorMsg = std::string("There is no MSC model called \"") + modelName + "\". Maybe it is not implemented or you misspelled the model name";
    Exception exc(errorMsg, __LINE__, __FILE__);
    exc.setFatal();
    throw exc;
  }
}


double MaterialEffects::effects(const std::vector< std::pair< MaterialProperties, M1x7 > >& points,
                                int materialsFXStart,
                                int materialsFXStop,
                                const double& mom,
                                const int& pdg,
                                M7x7* noise,
                                const M7x7* jacobian)
{

  if (materialInterface_ == nullptr) {
    std::string msg("MaterialEffects hasn't been initialized with a correct AbsMaterialInterface pointer!");
    std::runtime_error err(msg);
    throw err;
  }

  if (noEffects_) return 0.;

  bool doNoise(noise != nullptr);

  pdg_ = pdg;
  getParticleParameters(mom);

  double momLoss = 0.;

  for (auto it = points.begin() + materialsFXStart; it !=  points.begin() + materialsFXStop; ++it) { // loop over points

#ifdef DEBUG
    std::cerr << "     calculate matFX ";
    if (doNoise) std::cerr << " and noise";
    std::cerr << " for "; it->first.Print();
#endif

    double realPath = it->first.getSegmentLength();
    double stepSign(1.);
    if (realPath < 0)
      stepSign = -1.;
    realPath = fabs(realPath);

    if (realPath > 1.E-8) { // do material effects only if distance is not too small


      it->first.getMaterialProperties(matDensity_, matZ_, matA_, radiationLength_, mEE_);

      if (matZ_ > 1.E-3) { // don't calculate energy loss for vacuum

        if (energyLossBetheBloch_)
          momLoss += stepSign * this->energyLossBetheBloch(mom);
        if (doNoise && energyLossBetheBloch_ && noiseBetheBloch_)
          this->noiseBetheBloch(mom, *noise);

        if (doNoise && noiseCoulomb_)
          this->noiseCoulomb(mom, *noise, *((M1x3*) &it->second[3]) );

        if (energyLossBrems_)
          momLoss += stepSign * this->energyLossBrems(mom);
        if (doNoise && energyLossBrems_ && noiseBrems_)
          this->noiseBrems(mom, *noise);

      }
    }
  } // end loop over points

  return momLoss;
}


void MaterialEffects::stepper(const RKTrackRep* rep,
                              M1x7& state7,
                              const double& mom, // momentum
                              double& relMomLoss, // relative momloss for the step will be added
                              const int& pdg,
                              MaterialProperties& currentMaterial,
                              StepLimits& limits,
                              bool varField)
{

  static const double maxRelMomLoss = .005; // maximum relative momentum loss allowed
  static const double minStep = 1.E-4; // 1 µm

  // Trivial cases

  if (materialInterface_ == nullptr) {
    std::string msg("MaterialEffects hasn't been initialized with a correct AbsMaterialInterface pointer!");
    std::runtime_error err(msg);
    throw err;
  }

  if (noEffects_)
    return;

  if (relMomLoss > maxRelMomLoss) {
    limits.setLimit(stp_momLoss, 0);
    return;
  }


  double sMax = limits.getLowestLimitSignedVal(); // signed

  if (fabs(sMax) < minStep)
    return;



  pdg_ = pdg;
  getParticleParameters(mom);


  // make minStep
  state7[0] += limits.getStepSign() * minStep * state7[3];
  state7[1] += limits.getStepSign() * minStep * state7[4];
  state7[2] += limits.getStepSign() * minStep * state7[5];

  materialInterface_->initTrack(state7[0], state7[1], state7[2],
                                limits.getStepSign() * state7[3], limits.getStepSign() * state7[4], limits.getStepSign() * state7[5]);

  materialInterface_->getMaterialParameters(matDensity_, matZ_, matA_, radiationLength_, mEE_);
  currentMaterial.setMaterialProperties(matDensity_, matZ_, matA_, radiationLength_, mEE_);


#ifdef DEBUG
    std::cerr << "     currentMaterial "; currentMaterial.Print();
#endif

  // limit due to momloss
  double relMomLossPer_cm(0);
  stepSize_ = 1; // set stepsize for momLoss calculation

  if (matZ_ > 1.E-3) { // don't calculate energy loss for vacuum
    if (energyLossBetheBloch_) relMomLossPer_cm += this->energyLossBetheBloch(mom) / mom;
    if (energyLossBrems_)      relMomLossPer_cm += this->energyLossBrems(mom) / mom;
  }

  double maxStepMomLoss = (maxRelMomLoss - relMomLoss) / relMomLossPer_cm;
  limits.setLimit(stp_momLoss, maxStepMomLoss);

#ifdef DEBUG
    std::cerr << "     momLoss exceeded after a step of " <<  maxStepMomLoss << "\n";
#endif


  // now look for boundaries
  sMax = limits.getLowestLimitSignedVal();

  stepSize_ = limits.getStepSign() * minStep + materialInterface_->findNextBoundary(rep, state7, sMax, varField);
  if (fabs(stepSize_) < fabs(sMax)) {
    limits.setLimit(stp_boundary, stepSize_);
  }


  relMomLoss += relMomLossPer_cm * limits.getLowestLimitVal();
}


void MaterialEffects::getParticleParameters(double mom)
{
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(pdg_);
  charge_ = part->Charge() / (3.);
  mass_ = part->Mass();

  beta_ = mom / sqrt(mass_ * mass_ + mom * mom);

  //for numerical stability
  gammaSquare_ = 1. - beta_ * beta_;
  if (gammaSquare_ > 1.E-10) gammaSquare_ = 1. / gammaSquare_;
  else gammaSquare_ = 1.E10;
  gamma_ = sqrt(gammaSquare_);
}



//---- Energy-loss and Noise calculations -----------------------------------------

double MaterialEffects::energyLossBetheBloch(const double& mom)
{

  // calc dEdx_, also needed in noiseBetheBloch!
  dEdx_ = 0.307075 * matZ_ / matA_ * matDensity_ / (beta_ * beta_) * charge_ * charge_;
  double massRatio = me_ / mass_;
  double argument = gammaSquare_ * beta_ * beta_ * me_ * 1.E3 * 2. / ((1.E-6 * mEE_) * sqrt(1 + 2 * sqrt(gammaSquare_) * massRatio + massRatio * massRatio));
  if (argument <= exp(beta_ * beta_))
    dEdx_ = 0.;
  else {
    dEdx_ *= (log(argument) - beta_ * beta_); // Bethe-Bloch [MeV/cm]
    dEdx_ *= 1.E-3;  // in GeV/cm, hence 1.e-3
    if (dEdx_ < 0.) dEdx_ = 0;
  }

  double DE = fabs(stepSize_) * dEdx_; //always positive
  double momLoss = sqrt(mom * mom + 2.*sqrt(mom * mom + mass_ * mass_) * DE + DE * DE) - mom; //always positive

  //in vacuum it can numerically happen that momLoss becomes a small negative number.
  if (momLoss < 0.) return 0.;
  return momLoss;
}


void MaterialEffects::noiseBetheBloch(const double& mom,
                                      M7x7& noise) const
{

  // Code ported from GEANT 3

  // ENERGY LOSS FLUCTUATIONS; calculate sigma^2(E);
  double sigma2E = 0.;
  double zeta  = 153.4E3 * charge_ * charge_ / (beta_ * beta_) * matZ_ / matA_ * matDensity_ * fabs(stepSize_); // eV
  double Emax  = 2.E9 * me_ * beta_ * beta_ * gammaSquare_ / (1. + 2.*gamma_ * me_ / mass_ + (me_ / mass_) * (me_ / mass_)); // eV
  double kappa = zeta / Emax;

  if (kappa > 0.01) { // Vavilov-Gaussian regime
    sigma2E += zeta * Emax * (1. - beta_ * beta_ / 2.); // eV^2
  } else { // Urban/Landau approximation
    // calculate number of collisions Nc
    double I = 16. * pow(matZ_, 0.9); // eV
    double f2 = 0.;
    if (matZ_ > 2.) f2 = 2. / matZ_;
    double f1 = 1. - f2;
    double e2 = 10.*matZ_ * matZ_; // eV
    double e1 = pow((I / pow(e2, f2)), 1. / f1); // eV

    double mbbgg2 = 2.E9 * mass_ * beta_ * beta_ * gammaSquare_; // eV
    double Sigma1 = dEdx_ * 1.0E9 * f1 / e1 * (log(mbbgg2 / e1) - beta_ * beta_) / (log(mbbgg2 / I) - beta_ * beta_) * 0.6; // 1/cm
    double Sigma2 = dEdx_ * 1.0E9 * f2 / e2 * (log(mbbgg2 / e2) - beta_ * beta_) / (log(mbbgg2 / I) - beta_ * beta_) * 0.6; // 1/cm
    double Sigma3 = dEdx_ * 1.0E9 * Emax / (I * (Emax + I) * log((Emax + I) / I)) * 0.4; // 1/cm

    double Nc = (Sigma1 + Sigma2 + Sigma3) * fabs(stepSize_);

    if (Nc > 50.) { // truncated Landau distribution
      double sigmaalpha = 15.76;
      // calculate sigmaalpha  (see GEANT3 manual W5013)
      double RLAMED = -0.422784 - beta_ * beta_ - log(zeta / Emax);
      double RLAMAX =  0.60715 + 1.1934 * RLAMED + (0.67794 + 0.052382 * RLAMED) * exp(0.94753 + 0.74442 * RLAMED);
      // from lambda max to sigmaalpha=sigma (empirical polynomial)
      if (RLAMAX <= 1010.) {
        sigmaalpha =  1.975560
                      + 9.898841e-02 * RLAMAX
                      - 2.828670e-04 * RLAMAX * RLAMAX
                      + 5.345406e-07 * pow(RLAMAX, 3.)
                      - 4.942035e-10 * pow(RLAMAX, 4.)
                      + 1.729807e-13 * pow(RLAMAX, 5.);
      } else { sigmaalpha = 1.871887E+01 + 1.296254E-02 * RLAMAX; }
      // alpha=54.6  corresponds to a 0.9996 maximum cut
      if (sigmaalpha > 54.6) sigmaalpha = 54.6;
      sigma2E += sigmaalpha * sigmaalpha * zeta * zeta; // eV^2
    } else { // Urban model
      static const double alpha = 0.996;
      double Ealpha  = I / (1. - (alpha * Emax / (Emax + I))); // eV
      double meanE32 = I * (Emax + I) / Emax * (Ealpha - I); // eV^2
      sigma2E += fabs(stepSize_) * (Sigma1 * e1 * e1 + Sigma2 * e2 * e2 + Sigma3 * meanE32); // eV^2
    }
  }

  sigma2E *= 1.E-18; // eV -> GeV

  // update noise matrix
  noise[6 * 7 + 6] += (mom * mom + mass_ * mass_) / pow(mom, 6.) * sigma2E;
}


void MaterialEffects::noiseCoulomb(const double& mom,
                                   M7x7& noise,
                                   const M1x3& direction) const
{

//std::cerr << "MaterialEffects::noiseCoulomb" << std::endl;
//RKTools::printDim(jacobian,7,7);
  //double momTrue = 0.07;
  // MULTIPLE SCATTERING; calculate sigma^2
  double sigma2 = 0;
  assert(mscModelCode_ == 0 || mscModelCode_ == 1);
  const double step2 = stepSize_ * stepSize_;
  if (mscModelCode_ == 0) {// PANDA report PV/01-07 eq(43); linear in step length
    sigma2 = 225.E-6 * charge_ * charge_ / (beta_ * beta_ * mom * mom) * fabs(stepSize_) / radiationLength_ * matZ_ / (matZ_ + 1) * log(159.*pow(matZ_, -1. / 3.)) / log(287.*pow(matZ_, -0.5)); // sigma^2 = 225E-6*z^2/mom^2 * XX0/beta_^2 * Z/(Z+1) * ln(159*Z^(-1/3))/ln(287*Z^(-1/2)

  } else if (mscModelCode_ == 1) { //Highland not linear in step length formula taken from PDG book 2011 edition
    double stepOverRadLength = fabs(stepSize_) / radiationLength_;
    double logCor = (1 + 0.038 * log(stepOverRadLength));
    sigma2 = 0.0136 * 0.0136 * charge_ * charge_ / (beta_ * beta_ * mom * mom) * stepOverRadLength * logCor * logCor;
  }
  assert(sigma2 >= 0.0);
  //XXX std::cerr << "MaterialEffects::noiseCoulomb the MSC variance is " << sigma2 << std::endl;

  double noiseAfter[7 * 7]; // will hold the new MSC noise to cause by the current stepSize_ length
  memset(noiseAfter, 0x00, 7 * 7 * sizeof(double));

  double phi = atan2(direction[1], direction[0]);
  // cache sin and cos
  double sinTheta = sqrt(1 - direction[2] * direction[2]); // theta = arccos(direction[2])
  double cosTheta = direction[2];
  double sinPhi = sin(phi);
  double cosPhi = cos(phi);
  const double sinTheta2 = sinTheta * sinTheta;
  const double cosTheta2 = cosTheta * cosTheta;
  const double sinPhi2 = sinPhi * sinPhi;
  const double cosPhi2 = cosPhi * cosPhi;
  //this calculates the full projection of the MSC noise variance onto the 7D global coordinate system. Even taking into account the (co)variances of the position coordinates
  noiseAfter[0 * 7 + 0] =  sigma2 * step2 / 3.0 * (cosTheta2 + sinPhi2 * sinTheta2);
  noiseAfter[1 * 7 + 0] = -sigma2 * step2 / 3.0 * sinPhi * cosPhi * sinTheta2;
  noiseAfter[2 * 7 + 0] = -sigma2 * step2 / 3.0 * cosPhi * cosTheta * sinTheta;
  noiseAfter[3 * 7 + 0] =  sigma2 * stepSize_ * 0.5 * (cosTheta2 + sinPhi2 * sinTheta2);
  noiseAfter[4 * 7 + 0] = -sigma2 * stepSize_ * 0.5 * sinPhi * cosPhi * sinTheta2;
  noiseAfter[5 * 7 + 0] = -sigma2 * stepSize_ * 0.5 * cosPhi * cosTheta * sinTheta;
  noiseAfter[0 * 7 + 1] = noiseAfter[1 * 7 + 0];
  noiseAfter[1 * 7 + 1] = -sigma2 * step2 / 3.0 * (sinPhi2 * sinTheta2 - 1.0);
  noiseAfter[2 * 7 + 1] = -sigma2 * step2 / 3.0 * cosTheta * sinPhi * sinTheta;
  noiseAfter[3 * 7 + 1] = noiseAfter[4 * 7 + 0]; // Cov(x,a_y) = Cov(y,a_x)
  noiseAfter[4 * 7 + 1] = -sigma2 * stepSize_ * 0.5 * (sinPhi2 * sinTheta2 - 1.0);
  noiseAfter[5 * 7 + 1] = -sigma2 * stepSize_ * 0.5 * cosTheta * sinPhi * sinTheta;
  noiseAfter[0 * 7 + 2] = noiseAfter[2 * 7 + 0];
  noiseAfter[1 * 7 + 2] = noiseAfter[2 * 7 + 1];
  noiseAfter[2 * 7 + 2] =  sigma2 * step2 / 3.0 * sinTheta2;
  noiseAfter[3 * 7 + 2] = noiseAfter[5 * 7 + 0]; // Cov(z,a_x) = Cov(x,a_z)
  noiseAfter[4 * 7 + 2] = noiseAfter[5 * 7 + 1]; // Cov(y,a_z) = Cov(z,a_y)
  noiseAfter[5 * 7 + 2] =  sigma2 * stepSize_ * 0.5 * sinTheta2;
  noiseAfter[0 * 7 + 3] = noiseAfter[3 * 7 + 0];
  noiseAfter[1 * 7 + 3] = noiseAfter[3 * 7 + 1];
  noiseAfter[2 * 7 + 3] = noiseAfter[3 * 7 + 2];
  noiseAfter[3 * 7 + 3] =  sigma2 * (1.0 - sinTheta2 * cosPhi2);
  noiseAfter[4 * 7 + 3] = -sigma2 * sinTheta2 * cosPhi * sinPhi;
  noiseAfter[5 * 7 + 3] = -sigma2 * cosTheta * sinTheta * cosPhi;
  noiseAfter[0 * 7 + 4] = noiseAfter[4 * 7 + 0];
  noiseAfter[1 * 7 + 4] = noiseAfter[4 * 7 + 1];
  noiseAfter[2 * 7 + 4] = noiseAfter[4 * 7 + 2];
  noiseAfter[3 * 7 + 4] = noiseAfter[4 * 7 + 3];
  noiseAfter[4 * 7 + 4] =  sigma2 * (1.0 - sinTheta2 * sinPhi2);
  noiseAfter[5 * 7 + 4] = -sigma2 * cosTheta * sinTheta * sinPhi;
  noiseAfter[0 * 7 + 5] = noiseAfter[5 * 7 + 0];
  noiseAfter[1 * 7 + 5] = noiseAfter[5 * 7 + 1];
  noiseAfter[2 * 7 + 5] = noiseAfter[5 * 7 + 2];
  noiseAfter[3 * 7 + 5] = noiseAfter[5 * 7 + 3];
  noiseAfter[4 * 7 + 5] = noiseAfter[5 * 7 + 4];
  noiseAfter[5 * 7 + 5] = sigma2 * sinTheta2;
//    std::cerr << "new noise\n";
//    RKTools::printDim(noiseAfter, 7,7);
  for (unsigned int i = 0; i < 7 * 7; ++i) {
    noise[i] += noiseAfter[i];
  }
}


double MaterialEffects::energyLossBrems(const double& mom) const
{

  // Code ported from GEANT 3

  if (abs(pdg_) != 11) return 0; // only for electrons and positrons

#if !defined(BETHE)
  static const double C[101] = { 0.0, -0.960613E-01, 0.631029E-01, -0.142819E-01, 0.150437E-02, -0.733286E-04, 0.131404E-05, 0.859343E-01, -0.529023E-01, 0.131899E-01, -0.159201E-02, 0.926958E-04, -0.208439E-05, -0.684096E+01, 0.370364E+01, -0.786752E+00, 0.822670E-01, -0.424710E-02, 0.867980E-04, -0.200856E+01, 0.129573E+01, -0.306533E+00, 0.343682E-01, -0.185931E-02, 0.392432E-04, 0.127538E+01, -0.515705E+00, 0.820644E-01, -0.641997E-02, 0.245913E-03, -0.365789E-05, 0.115792E+00, -0.463143E-01, 0.725442E-02, -0.556266E-03, 0.208049E-04, -0.300895E-06, -0.271082E-01, 0.173949E-01, -0.452531E-02, 0.569405E-03, -0.344856E-04, 0.803964E-06, 0.419855E-02, -0.277188E-02, 0.737658E-03, -0.939463E-04, 0.569748E-05, -0.131737E-06, -0.318752E-03, 0.215144E-03, -0.579787E-04, 0.737972E-05, -0.441485E-06, 0.994726E-08, 0.938233E-05, -0.651642E-05, 0.177303E-05, -0.224680E-06, 0.132080E-07, -0.288593E-09, -0.245667E-03, 0.833406E-04, -0.129217E-04, 0.915099E-06, -0.247179E-07, 0.147696E-03, -0.498793E-04, 0.402375E-05, 0.989281E-07, -0.133378E-07, -0.737702E-02, 0.333057E-02, -0.553141E-03, 0.402464E-04, -0.107977E-05, -0.641533E-02, 0.290113E-02, -0.477641E-03, 0.342008E-04, -0.900582E-06, 0.574303E-05, 0.908521E-04, -0.256900E-04, 0.239921E-05, -0.741271E-07, -0.341260E-04, 0.971711E-05, -0.172031E-06, -0.119455E-06, 0.704166E-08, 0.341740E-05, -0.775867E-06, -0.653231E-07, 0.225605E-07, -0.114860E-08, -0.119391E-06, 0.194885E-07, 0.588959E-08, -0.127589E-08, 0.608247E-10};
  static const double xi = 2.51, beta = 0.99, vl = 0.00004;
#endif
#if defined(BETHE) // no MIGDAL corrections
  static const double C[101] = { 0.0, 0.834459E-02, 0.443979E-02, -0.101420E-02, 0.963240E-04, -0.409769E-05, 0.642589E-07, 0.464473E-02, -0.290378E-02, 0.547457E-03, -0.426949E-04, 0.137760E-05, -0.131050E-07, -0.547866E-02, 0.156218E-02, -0.167352E-03, 0.101026E-04, -0.427518E-06, 0.949555E-08, -0.406862E-02, 0.208317E-02, -0.374766E-03, 0.317610E-04, -0.130533E-05, 0.211051E-07, 0.158941E-02, -0.385362E-03, 0.315564E-04, -0.734968E-06, -0.230387E-07, 0.971174E-09, 0.467219E-03, -0.154047E-03, 0.202400E-04, -0.132438E-05, 0.431474E-07, -0.559750E-09, -0.220958E-02, 0.100698E-02, -0.596464E-04, -0.124653E-04, 0.142999E-05, -0.394378E-07, 0.477447E-03, -0.184952E-03, -0.152614E-04, 0.848418E-05, -0.736136E-06, 0.190192E-07, -0.552930E-04, 0.209858E-04, 0.290001E-05, -0.133254E-05, 0.116971E-06, -0.309716E-08, 0.212117E-05, -0.103884E-05, -0.110912E-06, 0.655143E-07, -0.613013E-08, 0.169207E-09, 0.301125E-04, -0.461920E-04, 0.871485E-05, -0.622331E-06, 0.151800E-07, -0.478023E-04, 0.247530E-04, -0.381763E-05, 0.232819E-06, -0.494487E-08, -0.336230E-04, 0.223822E-04, -0.384583E-05, 0.252867E-06, -0.572599E-08, 0.105335E-04, -0.567074E-06, -0.216564E-06, 0.237268E-07, -0.658131E-09, 0.282025E-05, -0.671965E-06, 0.565858E-07, -0.193843E-08, 0.211839E-10, 0.157544E-04, -0.304104E-05, -0.624410E-06, 0.120124E-06, -0.457445E-08, -0.188222E-05, -0.407118E-06, 0.375106E-06, -0.466881E-07, 0.158312E-08, 0.945037E-07, 0.564718E-07, -0.319231E-07, 0.371926E-08, -0.123111E-09};
  static const double xi = 2.10, beta_ = 1.00, vl = 0.001;
#endif

  double BCUT = 10000.; // energy up to which soft bremsstrahlung energy loss is calculated

  static const double THIGH = 100., CHIGH = 50.;
  double dedxBrems = 0.;

  if (BCUT > 0.) {
    double T, kc;

    if (BCUT >= mom) BCUT = mom; // confine BCUT to mom

    // T=mom,  confined to THIGH
    // kc=BCUT, confined to CHIGH ??
    if (mom >= THIGH) {
      T = THIGH;
      if (BCUT >= THIGH) kc = CHIGH;
      else kc = BCUT;
    } else {
      T = mom;
      kc = BCUT;
    }

    double E = T + me_; // total electron energy
    if (BCUT > T) kc = T;

    double X = log(T / me_);
    double Y = log(kc / (E * vl));

    double XX;
    int    K;
    double S = 0., YY = 1.;

    for (unsigned int I = 1; I <= 2; ++I) {
      XX = 1.;
      for (unsigned int J = 1; J <= 6; ++J) {
        K = 6 * I + J - 6;
        S = S + C[K] * XX * YY;
        XX = XX * X;
      }
      YY = YY * Y;
    }

    for (unsigned int I = 3; I <= 6; ++I) {
      XX = 1.;
      for (unsigned int J = 1; J <= 6; ++J) {
        K = 6 * I + J - 6;
        if (Y <= 0.) S = S + C[K] * XX * YY;
        else      S = S + C[K + 24] * XX * YY;
        XX = XX * X;
      }
      YY = YY * Y;
    }

    double SS = 0.;
    YY = 1.;

    for (unsigned int I = 1; I <= 2; ++I) {
      XX = 1.;
      for (unsigned int J = 1; J <= 5; ++J) {
        K = 5 * I + J + 55;
        SS = SS + C[K] * XX * YY;
        XX = XX * X;
      }
      YY = YY * Y;
    }

    for (unsigned int I = 3; I <= 5; ++I) {
      XX = 1.;
      for (unsigned int J = 1; J <= 5; ++J) {
        K = 5 * I + J + 55;
        if (Y <= 0.) SS = SS + C[K] * XX * YY;
        else      SS = SS + C[K + 15] * XX * YY;
        XX = XX * X;
      }
      YY = YY * Y;
    }

    S = S + matZ_ * SS;

    if (S > 0.) {
      double CORR = 1.;
#if !defined(BETHE)
      CORR = 1. / (1. + 0.805485E-10 * matDensity_ * matZ_ * E * E / (matA_ * kc * kc)); // MIGDAL correction factor
#endif

      double FAC = matZ_ * (matZ_ + xi) * E * E * pow((kc * CORR / T), beta) / (E + me_);
      if (FAC <= 0.) return 0.;
      dedxBrems = FAC * S;


      if (mom > THIGH) {
        double RAT;
        if (BCUT < THIGH) {
          RAT = BCUT / mom;
          S = (1. - 0.5 * RAT + 2.*RAT * RAT / 9.);
          RAT = BCUT / T;
          S = S / (1. - 0.5 * RAT + 2.*RAT * RAT / 9.);
        } else {
          RAT = BCUT / mom;
          S = BCUT * (1. - 0.5 * RAT + 2.*RAT * RAT / 9.);
          RAT = kc / T;
          S = S / (kc * (1. - 0.5 * RAT + 2.*RAT * RAT / 9.));
        }
        dedxBrems = dedxBrems * S; // GeV barn
      }

      dedxBrems = 0.60221367 * matDensity_ * dedxBrems / matA_; // energy loss dE/dx [GeV/cm]
    }
  }

  if (dedxBrems < 0.) dedxBrems = 0;

  double factor = 1.; // positron correction factor

  if (pdg_ == -11) {
    static const double AA = 7522100., A1 = 0.415, A3 = 0.0021, A5 = 0.00054;

    double ETA = 0.;
    if (matZ_ > 0.) {
      double X = log(AA * mom / matZ_ * matZ_);
      if (X > -8.) {
        if (X >= +9.) ETA = 1.;
        else {
          double W = A1 * X + A3 * pow(X, 3.) + A5 * pow(X, 5.);
          ETA = 0.5 + atan(W) / M_PI;
        }
      }
    }

    if (ETA < 0.0001) factor = 1.E-10;
    else if (ETA > 0.9999) factor = 1.;
    else {
      double E0 = BCUT / mom;
      if (E0 > 1.) E0 = 1.;
      if (E0 < 1.E-8) factor = 1.;
      else factor = ETA * (1. - pow(1. - E0, 1. / ETA)) / E0;
    }
  }

  double DE = fabs(stepSize_) * factor * dedxBrems; //always positive
  double momLoss = sqrt(mom * mom + 2.*sqrt(mom * mom + mass_ * mass_) * DE + DE * DE) - mom; //always positive

  return momLoss;
}


void MaterialEffects::noiseBrems(const double& mom,
                                 M7x7& noise) const
{

  // Code ported from GEANT 3

  if (abs(pdg_) != 11) return; // only for electrons and positrons

  double minusXOverLn2  = -1.442695 * fabs(stepSize_) / radiationLength_;
  double sigma2E = 1.44 * mom * mom * (pow(3., minusXOverLn2) - pow(4., minusXOverLn2));
  assert(sigma2E >= 0.0);
  noise[6 * 7 + 6] += (mom * mom + mass_ * mass_) / pow(mom, 6.) * sigma2E;

}

} /* End of namespace genfit */


