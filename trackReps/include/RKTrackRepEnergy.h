/* Copyright 2008-2015, Technische Universitaet Muenchen,
                        Ludwig-Maximilians-Universität München
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch
            & Tobias Schlüter

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

/** @addtogroup RKTrackRepEnergy
 * @{
 */

#ifndef genfit_RKTrackRepEnergy_h
#define genfit_RKTrackRepEnergy_h

#include "AbsTrackRep.h"
#include "StateOnPlane.h"
#include "RKTools.h"
#include "StepLimits.h"
#include "RKTrackRep.h" //RKStep

#include <algorithm>

namespace genfit {

/**
 * @brief AbsTrackRep with 5D track parameterization in plane coordinates: (q/p, u', v', u, v)
 *
 * q/p is charge over momentum.
 * u' and v' are direction tangents.
 * u and v are positions on a DetPlane.
 */
class RKTrackRepEnergy : public AbsTrackRep {


 public:

  RKTrackRepEnergy();
  RKTrackRepEnergy(int pdgCode, char propDir = 0);

  virtual ~RKTrackRepEnergy();

  virtual AbsTrackRep* clone() const {return new RKTrackRepEnergy(*this);}

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


  unsigned int getDim() const {return 5;}
  const char *getNameForLocalCoord(size_t i) {
    assert(0 <= i && i < getDim());
    switch (i) {
    case 0: return "qop";
    case 1: return "u'";
    case 2: return "v'";
    case 3: return "u";
    case 4: return "v";
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

  virtual double getRadiationLength() const;

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
              double& dlambda, M1x3& dT, RKMatrix<4,4>* pA) const;

  double RKstep(const M1x7& state7, const double h, const MaterialProperties& mat,
                M1x7& newState7, RKMatrix<7, 7>* pJ) const;

  //! The actual Runge Kutta propagation
  /** propagate state7 with step S. Fills SA (Start directions derivatives dA/S).
   *  This is a single Runge-Kutta step.
   *  If jacobian is NULL, only the state is propagated,
   *  otherwise also the 7x7 jacobian is calculated.
   *  If varField is false, the magnetic field will only be evaluated at the starting position.
   *  The return value is an estimation on how good the extrapolation is, and it is usually fine if it is > 1.
   *  It gives a suggestion how you must scale S so that the quality will be sufficient.
   */
  double RKPropagate(M1x7& state7,
                     M7x7* jacobian,
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

  void getState7(const StateOnPlane& state, M1x7& state7) const;
  void getState5(StateOnPlane& state, const SharedPlanePtr& plane, const M1x7& state7) const; // state7 must be on plane.

  void transformPM7(const MeasuredStateOnPlane& state,
                    M7x7& out7x7) const;

  void calcJ_pM_5x7(M5x7& J_pM, const TVector3& U, const TVector3& V, const M1x3& pTilde, double spu) const;

  void transformPM6(const MeasuredStateOnPlane& state,
                    M6x6& out6x6) const;

  void transformM7P(const M7x7& in7x7,
                    const M1x7& state7,
                    MeasuredStateOnPlane& state) const; // plane must already be set!

  void calcJ_Mp_7x5(M7x5& J_Mp, const TVector3& U, const TVector3& V, const M1x3& A) const;

  /**
   * Takes the 7x7 jacobian and noise for the transport of the global state startState7 
   * to the final global state destState7 and gives the 5x5 jacobian and noise, where 
   * the local coordinates are defined by startPlane and destPlane.
   */
  void projectJacobianAndNoise(const M1x7& startState7, const DetPlane& startPlane,
			       const M1x7& destState7, const DetPlane& destPlane,
			       const M7x7& jac, const M7x7& noise,
			       M5x5& jac5, M5x5& noise5) const;

  void transformM6P(const M6x6& in6x6,
                    const M1x7& state7,
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
              M1x7& state7,
              M7x7* jacobianT,
              double& coveredDistance, // signed
              double& flightTime,
              bool& checkJacProj,
              M7x7& noiseProjection,
              StepLimits& limits,
              bool onlyOneStep = false) const;

  double estimateStep(const M1x7& state7,
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
                M1x7& state7,
                double& flightTime,
                bool fillExtrapSteps,
                TMatrixDSym* cov = nullptr,
                bool onlyOneStep = false,
                bool stopAtBoundary = false,
                double maxStep = 1.E99) const;

  double momMag(const M1x7& state7) const;


  mutable StateOnPlane lastStartState_; //! state where the last extrapolation has started
  mutable StateOnPlane lastEndState_; //! state where the last extrapolation has ended
  mutable std::vector<RKStep> RKSteps_; //! RungeKutta steps made in the last extrapolation
  mutable int RKStepsFXStart_; //!
  mutable int RKStepsFXStop_; //!

  mutable TMatrixD fJacobian_; //!
  mutable TMatrixDSym fNoise_; //!

  void resetCache(const StateOnPlane& state) const;
  void checkCache(const StateOnPlane& state, const SharedPlanePtr& plane) const;
  mutable bool useCache_; //! use cached RKSteps_ for extrapolation
  mutable unsigned int cachePos_; //!

public:
  class propagator : public RKTrackRepEnergy::internalExtrapolator {
  public:
    propagator(const RKTrackRepEnergy* rep, const M1x7& state7, MaterialProperties& mat)
      : rep_(rep), state7_(state7), mat_(mat) { }
    void getInitialState(double posInitial[3], double dirInitial[3]) const {
      for(size_t i = 0; i < 3; ++i) {
        posInitial[i] = state7_[i];
        dirInitial[i] = state7_[i + 3];
      }
    }
    double extrapolateBy(double S, double posFinal[3], double dirFinal[3]) const {
      M1x3 SA;
      M1x7 state(state7_);
      double result = rep_->RKPropagate(state, 0, SA, S, mat_);
      for(size_t i = 0; i < 3; ++i) {
        posFinal[i] = state[i];
        dirFinal[i] = state[i + 3];
      }
      return result;
    }
    void moveStart(double posNew[3]) {
      std::copy(posNew, posNew + 3, state7_.vals);
    }
    void setMat(const MaterialProperties& mat) {
      mat_ = mat;
    }
  private:
    const RKTrackRepEnergy* rep_;
    M1x7 state7_;
    MaterialProperties mat_;
  };

 public:

  ClassDef(RKTrackRepEnergy, 1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTrackRepEnergy_h
