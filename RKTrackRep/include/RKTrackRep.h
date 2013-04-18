#ifndef genfit_RKTrackRep_h
#define genfit_RKTrackRep_h

#include <AbsTrackRep.h>


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


 private:

  mutable TVectorD lastStartState_; // state where the last extrapolation has started
  mutable TMatrixD jacobian_; // jacobian of the last extrapolation
  mutable TMatrixDSym noise_; // noise matrix of the last extrapolation
  mutable std::vector< MaterialProperties > materials_; // materials crossed in the last extrapolation

};

} /* End of namespace genfit */

#endif // genfit_RKTrackRep_h
