#ifndef genfit_AbsTrackRep_h
#define genfit_AbsTrackRep_h

#include <TVector3.h>

#include "MaterialInfo.h"
#include "MeasuredStateOnPlane.h"
#include "StateOnPlane.h"

namespace genfit {

  /** 
   *  Provides functionality to extrapolate a #GFStateOnPlane to another #GFDetPlane, or to the POCA to a line or a point.
   */
class AbsTrackRep {

 public:

  AbsTrackRep();
  AbsTrackRep(int pdgCode, char propDir = 0);

  virtual ~AbsTrackRep() {;}

  virtual AbsTrackRep* clone() const = 0;

  /** Extrapolates the stateInput to plane, and returns the extrapolation length
   * and, via reference, the extrapolated statePrediction.
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   */
  virtual double extrapolateToPlane(StateOnPlane* state,
      SharedPlanePtr plane,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToLine(StateOnPlane* state,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToPoint(StateOnPlane* state,
      const TVector3& point,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToCylinder(StateOnPlane* state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToSphere(StateOnPlane* state,
      double radius,
      const TVector3& point = TVector3(0.,0.,0.),
      bool stopAtBoundary = false) const = 0;

  /**
   * Use the Material information stored in the #TrackPoints
   */
  //virtual double extrapolateToTrackPoint() const;


  virtual unsigned int getDim() const = 0;

  virtual TVector3 getPos(const StateOnPlane* stateInput) const = 0;

  virtual TVector3 getMom(const StateOnPlane* stateInput) const = 0;
  virtual void getPosMom(const StateOnPlane* stateInput, TVector3& pos, TVector3& mom) const = 0;

  /** Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance */
  virtual void getPosMomCov(const MeasuredStateOnPlane* stateInput, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const = 0;

  int getPDG() const {return pdgCode_;}
  virtual double getCharge(const StateOnPlane* state) const = 0;
  char getPropDir() const {return propDir_;}

  /** Get the jacobian and noise matrix of the last extrapolation  */
  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise) const = 0;

  /** Get the jacobian and noise matrix of the last extrapolation if it would have been done in opposite direction  */
  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise) const = 0;


  virtual void setPosMom(StateOnPlane* stateInput, const TVector3& pos, const TVector3& mom) const = 0;
  virtual void setPosMomErr(MeasuredStateOnPlane* stateInput, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const = 0;
  virtual void setPosMomCov(MeasuredStateOnPlane* stateInput, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const = 0;

  //! Set propagation direction. (-1, 0, 1) -> (backward, auto, forward)
  void setPropDir(int dir) {
    if (dir>0) propDir_ = 1;
    else if (dir<0) propDir_ = -1;
    else propDir_ = 0;
  };

  //! Switch propagation direction. Has no effect if propDir_ is set to 0.
  void switchPropDir(){propDir_ = -1*propDir_;}


 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsTrackRep(const AbsTrackRep&) = default; // copy constructor
  AbsTrackRep& operator=(const AbsTrackRep&) = default; // assignment operator


  int pdgCode_;
  char propDir_;  // propagation direction (-1, 0, 1) -> (backward, auto, forward)

};

} /* End of namespace genfit */

#endif // genfit_AbsTrackRep_h
