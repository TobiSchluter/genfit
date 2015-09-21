/* Copyright 2008-2015, Technische Universitaet Muenchen,
                        Ludwig-Maximilians-Universität München
   Authors: Tobias Schlüter

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

/** @addtogroup RKTrackRepTime
 * @{
 */

#ifndef genfit_RKTrackRepTime_h
#define genfit_RKTrackRepTime_h

#include "AbsTrackRep.h"
#include "StateOnPlane.h"
#include "RKTools.h"
#include "StepLimits.h"
#include "RKTrackRep.h" //RKStep, ExtrapStep

#include <algorithm>

namespace genfit {

/**
 * @brief AbsTrackRep with 6D track parameterization in plane coordinates: (q/p, u', v', u, v, t)
 *
 * q/p is charge over momentum.
 * u' and v' are direction tangents.
 * u and v are positions on a DetPlane.
 * t is time with reference to some event time.
 */
class RKTrackRepTime : public AbsTrackRep {


 public:

  RKTrackRepTime();
  RKTrackRepTime(int pdgCode, char propDir = 0);

  virtual ~RKTrackRepTime();

  virtual AbsTrackRep* clone() const {return new RKTrackRepTime(*this);}

  virtual double extrapolateToPlane(StateOnPlane& state,
      const SharedPlanePtr& plane,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  using AbsTrackRep::extrapolateToLine;

  virtual double extrapolateToLine(StateOnPlane& state,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const {
    return extrapToPoint(state, point, NULL, stopAtBoundary, calcJacobianNoise);
  }

  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      const TMatrixDSym& G, // weight matrix (metric)
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const {
    return extrapToPoint(state, point, &G, stopAtBoundary, calcJacobianNoise);
  }

  virtual double extrapolateToCylinder(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  
  virtual double extrapolateToCone(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  virtual double extrapolateToSphere(StateOnPlane& state,
      double radius,
      const TVector3& point = TVector3(0.,0.,0.),
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  virtual double extrapolateBy(StateOnPlane& state,
      double step,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;


  unsigned int getDim() const {return 6;}
  const char *getNameForLocalCoord(size_t i) {
    assert(0 <= i && i < getDim());
    switch (i) {
    case 0: return "qop";
    case 1: return "u'";
    case 2: return "v'";
    case 3: return "u";
    case 4: return "v";
    case 5: return "t";
    }
    return 0; // Silence warnings about missing return statement.
  }


  virtual TVector3 getPos(const StateOnPlane& state) const;

  virtual TVector3 getMom(const StateOnPlane& state) const;
  virtual void getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const;

  virtual double getMomMag(const StateOnPlane& state) const;
  virtual double getMomVar(const MeasuredStateOnPlane& state) const;

  virtual TMatrixDSym get6DCov(const MeasuredStateOnPlane& state) const;
  virtual void getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const;
  virtual double getCharge(const StateOnPlane& state) const;
  virtual double getQop(const StateOnPlane& state) const {return state.getState()(0);}
  double getSpu(const StateOnPlane& state) const;
  double getTime(const StateOnPlane& state) const;

  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const;

  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const;

  std::vector<genfit::MatStep> getSteps() const;

  virtual double getRadiationLenght() const;

  virtual void setPosMom(StateOnPlane& state, const TVector3& pos, const TVector3& mom) const;
  virtual void setPosMom(StateOnPlane& state, const TVectorD& state6) const;
  virtual void setPosMomErr(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const;
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const;
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVectorD& state6, const TMatrixDSym& cov6x6) const;

  virtual void setChargeSign(StateOnPlane& state, double charge) const;
  virtual void setQop(StateOnPlane& state, double qop) const {state.getState()(0) = qop;}

  void setSpu(StateOnPlane& state, double spu) const;
  void setTime(StateOnPlane& state, double time) const;


  void derive(const double lambda, const M1x3& T,
              const double E, const double dEdx, const double d2EdxdE, const double B[3],
              double& dlambda, M1x3& dT, double& dTime, M5x5* pA) const;

  double RKstep(const M1x8& stateGlobal, const double h, const MaterialProperties& mat,
                M1x8& newStateGlobal, M8x8* pJ) const;

  //! The actual Runge Kutta propagation
  /** propagate state7 with step S. Fills SA (Start directions derivatives dA/S).
   *  This is a single Runge-Kutta step.
   *  If jacobian is NULL, only the state is propagated,
   *  otherwise also the 7x7 jacobian is calculated.
   *  If varField is false, the magnetic field will only be evaluated at the starting position.
   *  The return value is an estimation on how good the extrapolation is, and it is usually fine if it is > 1.
   *  It gives a suggestion how you must scale S so that the quality will be sufficient.
   */
  double RKPropagate(M1x8& stateGlobal,
                     M8x8* jacobian,
                     M1x3& SA,
                     double S,
                     const MaterialProperties& mat) const;

  virtual bool isSameType(const AbsTrackRep* other);
  virtual bool isSame(const AbsTrackRep* other);

 private:

  void initArrays() const;

  virtual double extrapToPoint(StateOnPlane& state,
      const TVector3& point,
      const TMatrixDSym* G = NULL, // weight matrix (metric)
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const;

  void getStateGlobal(const StateOnPlane& stateLocal, M1x8& stateGlobal) const;
  void getStateLocal(StateOnPlane& stateLocal, const M1x8& stateGlobal) const; // state8 must already lie on plane of state!

  void transformPM8(const MeasuredStateOnPlane& state,
                    M8x8& out8x8) const;

  void calcJ_pM_6x8(M6x8& J_pM, const TVector3& U, const TVector3& V, const M1x3& pTilde, double spu) const;

  void transformPM6(const MeasuredStateOnPlane& state,
                    M6x6& out6x6) const;

  void transformM8P(const M8x8& in8x8,
                    const M1x8& state8,
                    MeasuredStateOnPlane& state) const; // plane must already be set!

  void calcJ_Mp_8x6(M8x6& J_Mp, const TVector3& U, const TVector3& V, const M1x3& A) const;

  void calcForwardJacobianAndNoise(const M1x8& startState8, const DetPlane& startPlane,
				   const M1x8& destState8, const DetPlane& destPlane) const;

  void transformM6P(const M6x6& in6x6,
                    const M1x8& stateGlobal,
                    MeasuredStateOnPlane& state) const; // plane and charge must already be set!

  //! Propagates the particle through the magnetic field.
  /** If the propagation is successful and the plane is reached, the function returns true.
    * Propagated state and the jacobian of the extrapolation are written to state7 and jacobianT.
    * The jacobian is only calculated if jacobianT != NULL.
    * In the main loop of the Runge Kutta algorithm, the estimateStep() is called
    * and may reduce the estimated stepsize so that a maximum momentum loss will not be exceeded,
    * and stop at material boundaries.
    * If this is the case, RKutta() will only propagate the reduced distance and then return. This is to ensure that
    * material effects, which are calculated after the propagation, are taken into account properly.
    */
  bool RKutta(const M1x4& SU,
              const DetPlane& plane,
              double charge,
              double mass,
              M1x8& stateGlobal,
              M8x8* jacobianT,
              double& coveredDistance, // signed
              double& flightTime,
              bool& checkJacProj,
              M7x7& noiseProjection,
              StepLimits& limits,
              bool onlyOneStep = false) const;

  double estimateStep(const M1x8& stateGlobal,
                      const M1x4& SU,
                      const DetPlane& plane,
                      const double& charge,
                      double& relMomLoss,
                      StepLimits& limits,
                      MaterialProperties& matForStep) const;

  TVector3 pocaOnLine(const TVector3& linePoint,
                     const TVector3& lineDirection,
                     const TVector3& point) const;

  //! Handles propagation and material effects
  /** #extrapolateToPlane(), #extrapolateToPoint() and #extrapolateToLine() etc. call this function.
    * #Extrap() needs a plane as an argument, hence #extrapolateToPoint() and #extrapolateToLine() create virtual detector planes.
    * In this function, #RKutta() is called and the resulting points and point paths are filtered
    * so that the direction doesn't change and tiny steps are filtered out.
    * After the propagation the material effects are called via the MaterialEffects singleton.
    * #Extrap() will loop until the plane is reached, unless the propagation fails or the maximum number of
    * iterations is exceeded.
    */
  double Extrap(const DetPlane& startPlane, // plane where Extrap starts
                const DetPlane& destPlane, // plane where Extrap has to extrapolate to
                double charge,
                double mass,
                bool& isAtBoundary,
                M1x8& state8,
                double& flightTime,
                bool fillExtrapSteps,
                TMatrixDSym* cov = nullptr,
                bool onlyOneStep = false,
                bool stopAtBoundary = false,
                double maxStep = 1.E99) const;

  void checkCache(const StateOnPlane& state, const SharedPlanePtr* plane) const;

  double momMag(const M1x8& stateGlobal) const;


  mutable StateOnPlane lastStartState_; //! state where the last extrapolation has started
  mutable StateOnPlane lastEndState_; //! state where the last extrapolation has ended
  mutable std::vector<TRKStep<8> > RKSteps_; //! RungeKutta steps made in the last extrapolation
  mutable int RKStepsFXStart_; //!
  mutable int RKStepsFXStop_; //!
  mutable std::vector<TExtrapStep<8> > ExtrapSteps_; //! steps made in Extrap during last extrapolation

  mutable TMatrixD fJacobian_; //!
  mutable TMatrixDSym fNoise_; //!

  mutable bool useCache_; //! use cached RKSteps_ for extrapolation
  mutable unsigned int cachePos_; //!

  // auxiliary variables and arrays
  // needed in Extrap()
  mutable StepLimits limits_; //!
  mutable M7x7 noiseProjection_; //!
public:
  class propagator : public RKTrackRepTime::internalExtrapolator {
  public:
    propagator(const RKTrackRepTime* rep, const M1x8& stateGlobal, MaterialProperties& mat)
      : rep_(rep), stateGlobal_(stateGlobal), mat_(mat) {}
    void getInitialState(double posInitial[3], double dirInitial[3]) const {
      for(size_t i = 0; i < 3; ++i) {
        posInitial[i] = stateGlobal_[i];
        dirInitial[i] = stateGlobal_[i + 3];
      }
    }
    double extrapolateBy(double S, double posFinal[3], double dirFinal[3]) const {
      M1x3 SA;
      M1x8 state;
      state = stateGlobal_;
      double result = rep_->RKPropagate(state, 0, SA, S, mat_);
      for(size_t i = 0; i < 3; ++i) {
        posFinal[i] = state[i];
        dirFinal[i] = state[i + 3];
      }
      return result;
    }
    void moveStart(double posNew[3]) {
      std::copy(posNew, posNew + 3, stateGlobal_.vals);
    }
    void setMat(const MaterialProperties& mat) {
      mat_ = mat;
    }
  private:
    const RKTrackRepTime* rep_;
    M1x8 stateGlobal_;
    MaterialProperties mat_;
  };

 public:

  ClassDef(RKTrackRepTime, 1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTrackRepTime_h
