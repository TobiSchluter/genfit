#ifndef genfit_AbsTrackRep_h
#define genfit_AbsTrackRep_h

#include <TVector3.h>
#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include "SharedPlanePtr.h"
//#include "MaterialInfo.h"

namespace genfit {

class StateOnPlane;
class MeasuredStateOnPlane;

  /** 
   *  Provides functionality to extrapolate a #GFStateOnPlane to another #GFDetPlane, or to the POCA to a line or a point.
   */
class AbsTrackRep : public TObject {

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
  TVector3 getDir(const StateOnPlane* stateInput) const {return getMom(stateInput).Unit();}
  virtual void getPosMom(const StateOnPlane* stateInput, TVector3& pos, TVector3& mom) const = 0;
  void getPosDir(const StateOnPlane* stateInput, TVector3& pos, TVector3& dir) const {getPosMom(stateInput, pos, dir); dir.SetMag(1.);}
  virtual TVectorD get6DState(const StateOnPlane* stateInput) const;

  /** Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance */
  virtual void getPosMomCov(const MeasuredStateOnPlane* stateInput, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const = 0;
  virtual void get6DStateCov(const MeasuredStateOnPlane* stateInput, TVectorD& stateVec, TMatrixDSym& cov) const;

  int getPDG() const {return pdgCode_;}
  virtual double getCharge(const StateOnPlane* state) const = 0;
  char getPropDir() const {return propDir_;}

  /** Get the jacobian and noise matrix of the last extrapolation  */
  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise) const = 0;

  /** Get the jacobian and noise matrix of the last extrapolation if it would have been done in opposite direction  */
  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise) const = 0;

  /** Calculate Jacobian of transportation numerically. Slow but accurate. Can be used to validate (semi)analytic calculations. */
  void calcJacobianNumerically(const genfit::StateOnPlane* origState,
                                   const genfit::SharedPlanePtr destPlane,
                                   TMatrixD& jacobian);

  virtual void setPosMom(StateOnPlane* stateInput, const TVector3& pos, const TVector3& mom) const = 0;
  virtual void setPosMom(StateOnPlane* stateInput, const TVectorD& state6) const = 0;
  virtual void setPosMomErr(MeasuredStateOnPlane* stateInput, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const = 0;
  virtual void setPosMomCov(MeasuredStateOnPlane* stateInput, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const = 0;
  virtual void setPosMomCov(MeasuredStateOnPlane* stateInput, const TVectorD& state6, const TMatrixDSym& cov6x6) const = 0;

  //! Set propagation direction. (-1, 0, 1) -> (backward, auto, forward)
  void setPropDir(int dir) {
    if (dir>0) propDir_ = 1;
    else if (dir<0) propDir_ = -1;
    else propDir_ = 0;
  };

  //! Switch propagation direction. Has no effect if propDir_ is set to 0.
  void switchPropDir(){propDir_ = -1*propDir_;}

  virtual void Print(const Option_t* = "") const;

 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsTrackRep(const AbsTrackRep&); // copy constructor
  AbsTrackRep& operator=(const AbsTrackRep&); // assignment operator


  int pdgCode_;
  char propDir_;  // propagation direction (-1, 0, 1) -> (backward, auto, forward)

};

} /* End of namespace genfit */

#endif // genfit_AbsTrackRep_h
