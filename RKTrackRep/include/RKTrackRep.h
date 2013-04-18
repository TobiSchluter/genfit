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

  virtual AbsTrackRep* clone() const = 0;

  /** Extrapolates the stateInput to plane, and returns the extrapolation length
   * and, via reference, the extrapolated statePrediction.
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   */
  virtual double extrapolateToPlane(const StateOnPlane& stateInput,
      StateOnPlane& statePrediction,
      sharedPlanePtr plane,
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

  virtual double getCharge() const override;

  /** Get the jacobian of the last extrapolation  */
  virtual TMatrixD getForwardJacobian() const override;

  /** Get the jacobian of the last extrapolation if it would have been done in opposite direction  */
  virtual TMatrixD getBackwardJacobian() const override;

  /** Get the noise matrix of the last extrapolation  */
  virtual TMatrixDSym getForwardNoise() const override;

  /** Get the noise matrix of the last extrapolation if it would have been done in opposite direction  */
  virtual TMatrixDSym getBackwardNoise() const override;


  virtual void setPosMom(StateOnPlane* stateInput, const TVector3& pos, const TVector3& mom) const override;
  virtual void setPosMomCov(MeasuredStateOnPlane* stateInput, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov) const override;


  //! The actual Runge Kutta propagation
  /** propagate #state7 with step #S. Fills #SA (Start directions derivatives dA/S).
   *  If #cov is NULL, only the state is propagated,
   *  otherwise also the 7x7 jacobian (#cov) is calculated.
   *  If #varField is false, the magnetic field will only be evaluated at the starting position.
   *  The return value is an estimation on how good the extrapolation is, and it is usually fine if it is > 1.
   *  It gives a suggestion how you must scale #S so that the quality will be sufficient.
   */
  double RKPropagate(M1x7& state7,
                     M7x7* cov,
                     M1x3& SA,
                     double S,
                     bool varField = true) const;


 private:

  void calcStateCov(const TVector3& pos,
                    const TVector3& mom,
                    const TVector3& poserr,
                    const TVector3& momerr);

  void calcState(const TVector3& pos,
                 const TVector3& mom);

  void getState7(M1x7& state7);
  void getState7(M1x7& state7, const TVectorD& state5, const DetPlane& pl, const double& spu);
  TVectorD getState5(const M1x7& state7, const DetPlane& pl, double& spu);

  void transformPM7(const TMatrixD& in5x5,
                    M7x7& out7x7,
                    const DetPlane& pl,
                    const TVectorD& state5,
                    const double& spu,
                    TMatrixD* Jac = NULL);

  void transformPM6(const TMatrixDSym& in5x5,
                    M6x6& out6x6,
                    const DetPlane& pl,
                    const TVectorD& state5,
                    const double& spu,
                    TMatrixD* Jac = NULL);

  void transformM7P(const M7x7& in7x7,
                    TMatrixDSym& out5x5,
                    const DetPlane& pl,
                    const M1x7& state7,
                    TMatrixD* Jac = NULL);

  void transformM6P(const M6x6& in6x6,
                    TMatrixDSym& out5x5,
                    const DetPlane& pl,
                    const M1x7& state7,
                    TMatrixD* Jac = NULL);

  //! Propagates the particle through the magnetic field.
  /** If the propagation is successfull and the plane is reached, the function returns true.
    * Propagated state and the jacobian of the extrapolation are written to #state7 and #cov.
    * The jacobian is only calculated if #cov != NULL.
    * In the main loop of the Runge Kutta algorithm, the #estimateStep() is called
    * and may reduce the estimated stepsize so that a maximum momentum loss will not be exceeded.
    * If this is the case, #RKutta() will only propagate the reduced distance and then return. This is to ensure that
    * material effects, which are calculated after the propagation, are taken into account properly.
    */
  bool RKutta (const DetPlane& plane,
               M1x7& state7,
               M7x7* cov,
               double& coveredDistance,
               std::vector<MaterialProperties>& points,
               bool& checkJacProj,
               TMatrixD& noiseProjection,
               bool onlyOneStep = false,
               double maxStep = 1.E99);

  double estimateStep(std::vector<MaterialProperties>& points,
                      const TVector3& pos,
                      const TVector3& dir,
                      const M1x4& SU,
                      const DetPlane& plane,
                      const double& mom,
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
                M1x7& state7,
                M7x7* cov=NULL,
                bool onlyOneStep = false,
                double maxStep = 1.E99);



  mutable TVectorD lastStartState_; // state where the last extrapolation has started
  mutable TMatrixD jacobian_; // jacobian of the last extrapolation
  mutable TMatrixDSym noise_; // noise matrix of the last extrapolation
  mutable std::vector< MaterialProperties > materials_; // materials crossed in the last extrapolation

};

} /* End of namespace genfit */

#endif // genfit_RKTrackRep_h
