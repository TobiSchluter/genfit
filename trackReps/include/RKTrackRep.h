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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_RKTrackRep_h
#define genfit_RKTrackRep_h

#include "AbsTrackRep.h"
#include "StateOnPlane.h"
#include "MaterialInfo.h"

#include "RKTools.h"
#include "StepLimits.h"


namespace genfit {

struct RKStep {
  MaterialProperties materialProperties_;
  M1x7 state7_; // 7D state vector
};

struct ExtrapStep {
  M5x5 jac_; // 5D jacobian of transport
  M5x5 noise_; // 5D noise matrix
};

class RKTrackRep : public AbsTrackRep {

  /**
   * state5: (q/p, u', v'. u. v)
   */

 public:

  RKTrackRep();
  RKTrackRep(int pdgCode, char propDir = 0);

  virtual ~RKTrackRep();

  virtual AbsTrackRep* clone() const {return new RKTrackRep(*this);}

  /** Extrapolates the stateInput to plane, and returns the extrapolation length
   * and, via reference, the extrapolated statePrediction.
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   */
  virtual double extrapolateToPlane(StateOnPlane* state,
      SharedPlanePtr plane,
      bool stopAtBoundary = false) const;

  virtual double extrapolateToLine(StateOnPlane* state,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false) const;

  virtual double extrapolateToPoint(StateOnPlane* state,
      const TVector3& point,
      bool stopAtBoundary = false) const;

  virtual double extrapolateToCylinder(StateOnPlane* state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false) const;

  virtual double extrapolateToSphere(StateOnPlane* state,
      double radius,
      const TVector3& point = TVector3(0.,0.,0.),
      bool stopAtBoundary = false) const;

  /**
   * Use the Material information stored in the #TrackPoints
   */
  //virtual double extrapolateToTrackPoint() const;


  unsigned int getDim() const {return 5;}

  virtual TVector3 getPos(const StateOnPlane* stateInput) const;

  virtual TVector3 getMom(const StateOnPlane* stateInput) const;
  virtual void getPosMom(const StateOnPlane* stateInput, TVector3& pos, TVector3& mom) const;

  /** Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance */
  virtual void getPosMomCov(const MeasuredStateOnPlane* stateInput, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const;
  virtual double getCharge(const StateOnPlane* state) const;
  double getSpu(const StateOnPlane* state) const;

  /** Get the jacobian and noise matrix of the last extrapolation  */
  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise) const;

  /** Get the jacobian and noise matrix of the last extrapolation if it would have been done in opposite direction  */
  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise) const;


  virtual void setPosMom(StateOnPlane* state, const TVector3& pos, const TVector3& mom) const;
  virtual void setPosMom(StateOnPlane* stateInput, const TVectorD& state6) const;
  virtual void setPosMomErr(MeasuredStateOnPlane* state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const;
  virtual void setPosMomCov(MeasuredStateOnPlane* state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const;
  virtual void setPosMomCov(MeasuredStateOnPlane* state, const TVectorD& state6, const TMatrixDSym& cov6x6) const;


  void setCharge(StateOnPlane* state, double charge) const {(state->getAuxInfo())(0) = charge;}
  void setSpu(StateOnPlane* state, double spu) const {(state->getAuxInfo())(1) = spu;}

  //! The actual Runge Kutta propagation
  /** propagate #state7 with step #S. Fills #SA (Start directions derivatives dA/S).
   *  If #cov is NULL, only the state is propagated,
   *  otherwise also the 7x7 jacobian (#cov) is calculated.
   *  If #varField is false, the magnetic field will only be evaluated at the starting position.
   *  The return value is an estimation on how good the extrapolation is, and it is usually fine if it is > 1.
   *  It gives a suggestion how you must scale #S so that the quality will be sufficient.
   */
  double RKPropagate(M1x7& state7,
                     M7x7* jacobian,
                     M1x3& SA,
                     double S,
                     bool varField = true,
                     bool calcOnlyLastRowOfJ = false) const;


 private:

  void initArrays() const;

  void getState7(const StateOnPlane* state, M1x7& state7) const;
  void getState5(StateOnPlane* state, const M1x7& state7) const; // state7 must already lie on plane of state!

  void transformPM7(const MeasuredStateOnPlane* state,
                    M7x7& out7x7) const;

  void calcJ_pM_5x7(const TVector3& U, const TVector3& V, const M1x3& pTilde, double spu) const;

  void transformPM6(const MeasuredStateOnPlane* state,
                    M6x6& out6x6) const;

  void transformM7P(const M7x7& in7x7,
                    const M1x7& state7,
                    MeasuredStateOnPlane* state) const; // plane must already be set!

  void calcJ_Mp_7x5(const TVector3& U, const TVector3& V, const TVector3& W, const M1x3& A) const;

  void transformM6P(const M6x6& in6x6,
                    const M1x7& state7,
                    MeasuredStateOnPlane* state) const; // plane and charge must already be set!

  //! Propagates the particle through the magnetic field.
  /** If the propagation is successful and the plane is reached, the function returns true.
    * Propagated state and the jacobian of the extrapolation are written to #state7 and #jacobian.
    * The jacobian is only calculated if #jacobian != NULL.
    * In the main loop of the Runge Kutta algorithm, the #estimateStep() is called
    * and may reduce the estimated stepsize so that a maximum momentum loss will not be exceeded.
    * If this is the case, #RKutta() will only propagate the reduced distance and then return. This is to ensure that
    * material effects, which are calculated after the propagation, are taken into account properly.
    */
  bool RKutta(const M1x4& SU,
              const DetPlane& plane,
              double charge,
              M1x7& state7,
              M7x7* jacobianT,
              double& coveredDistance,
              bool& checkJacProj,
              TMatrixD& noiseProjection,
              StepLimits& limits,
              bool onlyOneStep = false,
              bool calcOnlyLastRowOfJ = false) const;

  double estimateStep(const M1x7& state7,
                      const M1x4& SU,
                      const DetPlane& plane,
                      const double& charge,
                      double& relMomLoss,
                      StepLimits& limits) const;

  TVector3 pocaOnLine(const TVector3& linePoint,
                     const TVector3& lineDirection,
                     const TVector3& point) const;

  //! Handles propagation and material effects
  /** #extrapolate(), #extrapolateToPoint() and #extrapolateToLine() call this function.
    * #Extrap() needs a plane as an argument, hence #extrapolateToPoint() and #extrapolateToLine() create virtual detector planes.
    * In this function, #RKutta() is called and the resulting points and point paths are filtered
    * so that the direction doesn't change and tiny steps are filtered out.
    * After the propagation the material effects are called via the #MaterialEffects singleton.
    * #Extrap() will loop until the plane is reached, unless the propagation fails or the maximum number of
    * iterations is exceeded.
    * #fXX0 is also updated here.
    */
  double Extrap(const DetPlane& startPlane, // plane where Extrap starts
                const DetPlane& destPlane, // plane where Extrap has to extrapolate to
                double charge,
                bool& isAtBoundary,
                M1x7& state7,
                TMatrixDSym* cov = nullptr,
                bool onlyOneStep = false,
                bool stopAtBoundary = false,
                double maxStep = 1.E99) const;

  void checkCache(const StateOnPlane* state) const;



  mutable StateOnPlane lastStartState_; //! state where the last extrapolation has started
  mutable std::vector<RKStep> RKSteps_; //! RungeKutta steps made in the last extrapolation
  mutable int RKStepsFXStart_; //!
  mutable int RKStepsFXStop_; //!
  mutable std::vector<ExtrapStep> ExtrapSteps_; //! steps made in Extrap during last extrapolation

  mutable bool useCache_; //! use cached RKSteps_ for extrapolation

  // auxiliary variables and arrays
  // needed in Extrap()
  mutable M7x7 noiseArray_; //! noise matrix of the last extrapolation
  mutable M7x7 J_MMT_; //!
  // needed in transform...
  mutable M5x7 J_pM_5x7_; //!  // FIXME this is actually (J_Mp)^T
  mutable M5x6 J_pM_5x6_; //!  // FIXME this is actually (J_Mp)^T
  mutable M7x5 J_Mp_7x5_; //!  // FIXME this is actually (J_pM)^T
  mutable M6x5 J_Mp_6x5_; //!  // FIXME this is actually (J_pM)^T

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTrackRep_h
