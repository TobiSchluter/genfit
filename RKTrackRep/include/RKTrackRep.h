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

/** @addtogroup genfit
 * @{
 */

#ifndef genfit_RKTrackRep_h
#define genfit_RKTrackRep_h

#include <AbsTrackRep.h>

#include "RKTools.h"


namespace genfit {

class RKTrackRep : public AbsTrackRep {

 public:

  RKTrackRep();
  RKTrackRep(int pdgCode, char propDir = 0);

  virtual ~RKTrackRep() {;}

  virtual AbsTrackRep* clone() const {return new RKTrackRep(*this);}

  /** Extrapolates the stateInput to plane, and returns the extrapolation length
   * and, via reference, the extrapolated statePrediction.
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   */
  virtual double extrapolateToPlane(const StateOnPlane& stateInput,
      StateOnPlane& statePrediction,
      SharedPlanePtr plane,
      bool stopAtBoundary = false) const override;

  virtual double extrapolateToLine(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false) const override;

  virtual double extrapolateToPoint(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& point,
      bool stopAtBoundary = false) const override;

  virtual double extrapolateToCylinder(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      double radius,
      bool stopAtBoundary = false) const override;

  virtual double extrapolateToSphere(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& point,
      double radius,
      bool stopAtBoundary = false) const override;

  /**
   * Use the Material information stored in the #TrackPoints
   */
  //virtual double extrapolateToTrackPoint() const;


  virtual TVector3 getPos(const StateOnPlane* stateInput) const override;

  virtual TVector3 getMom(const StateOnPlane* stateInput) const override;
  virtual void getPosMom(const StateOnPlane* stateInput, TVector3& pos, TVector3& mom) const override;

  /** Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance */
  virtual void getPosMomCov(const MeasuredStateOnPlane* stateInput, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const override;

  virtual double getCharge(const StateOnPlane* state) const override {return (state->getAuxInfo())(0);}
  double getSpu(const StateOnPlane* state) const {return (state->getAuxInfo())(1);}

  /** Get the jacobian of the last extrapolation  */
  virtual TMatrixD getForwardJacobian() const override;

  /** Get the jacobian of the last extrapolation if it would have been done in opposite direction  */
  virtual TMatrixD getBackwardJacobian() const override;

  /** Get the noise matrix of the last extrapolation  */
  virtual TMatrixDSym getForwardNoise() const override;

  /** Get the noise matrix of the last extrapolation if it would have been done in opposite direction  */
  virtual TMatrixDSym getBackwardNoise() const override;


  virtual void setPosMom(StateOnPlane* state, const TVector3& pos, const TVector3& mom) const override;
  virtual void setPosMomErr(MeasuredStateOnPlane* state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const override;
  virtual void setPosMomCov(MeasuredStateOnPlane* state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const override;


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
                     bool varField = true) const;


 private:

  void initArrays();

  void getState7(const StateOnPlane* state, M1x7& state7) const;
  void getState5(StateOnPlane* state, const M1x7& state7) const;

  void transformPM7(const MeasuredStateOnPlane* state,
                    M7x7& out7x7,
                    TMatrixD* Jac = NULL) const;

  void transformPM6(const MeasuredStateOnPlane* state,
                    M6x6& out6x6,
                    TMatrixD* Jac = NULL) const;

  void transformM7P(const M7x7& in7x7,
                    const M1x7& state7,
                    MeasuredStateOnPlane* state, // plane must already be set!
                    TMatrixD* Jac = NULL) const;

  void transformM6P(const M6x6& in6x6,
                    const M1x7& state7,
                    MeasuredStateOnPlane* state, // plane and charge must already be set!
                    TMatrixD* Jac = NULL) const;

  //! Propagates the particle through the magnetic field.
  /** If the propagation is successful and the plane is reached, the function returns true.
    * Propagated state and the jacobian of the extrapolation are written to #state7 and #jacobian.
    * The jacobian is only calculated if #jacobian != NULL.
    * In the main loop of the Runge Kutta algorithm, the #estimateStep() is called
    * and may reduce the estimated stepsize so that a maximum momentum loss will not be exceeded.
    * If this is the case, #RKutta() will only propagate the reduced distance and then return. This is to ensure that
    * material effects, which are calculated after the propagation, are taken into account properly.
    */
  bool RKutta(const DetPlane& plane,
              double charge,
              M1x7& state7,
              M7x7* jacobian,
              double& coveredDistance,
              bool& checkJacProj,
              TMatrixD& noiseProjection,
              bool onlyOneStep = false,
              double maxStep = 1.E99) const;

  double estimateStep(const M1x7& state7,
                      const M1x4& SU,
                      const DetPlane& plane,
                      const double& charge,
                      double& relMomLoss,
                      bool& momLossExceeded,
                      bool& atPlane,
                      double maxStep = 1.E99) const;

  TVector3 poca2Line(const TVector3& extr1,
                     const TVector3& extr2,
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
  double Extrap(const DetPlane& plane,
                double charge,
                M1x7& state7,
                M7x7* cov=NULL,
                bool onlyOneStep = false,
                double maxStep = 1.E99) const;



  mutable StateOnPlane lastStartState_; //! state where the last extrapolation has started
  mutable TMatrixD jacobian_; //! jacobian of the last extrapolation
  mutable TMatrixDSym noise_; //! noise matrix of the last extrapolation
  mutable std::vector< std::pair< MaterialProperties, M1x7 > > materials_; //! materials crossed in the last extrapolation, together with 7D states at start of each step


  // auxiliary variables and arrays
  // needed in Extrap()
  mutable M7x7 fNoise; //!
  mutable M7x7 fOldCov; //!
  // needed in transform...
  mutable M5x7 fJ_pM_5x7; //!
  mutable M5x6 fJ_pM_5x6; //!
  mutable M7x5 fJ_Mp_7x5; //!
  mutable M6x5 fJ_Mp_6x5; //!

};

} /* End of namespace genfit */

#endif // genfit_RKTrackRep_h
