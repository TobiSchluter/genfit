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

  /** Extrapolates the state to plane, and returns the extrapolation length
   * and, via reference, the extrapolated statePrediction.
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   * If state has a covariance, jacobian and noise matrices will be calculated and the covariance will be propagated.
   * If state has no covariance, jacobian and noise will only be calculated if calcJacobianNoise == true.
   */
  virtual double extrapolateToPlane(StateOnPlane& state,
      const SharedPlanePtr& plane,
      bool stopAtBoundary = false,
      bool calcJacobianNoise = false) const = 0;

  virtual double extrapolateToLine(StateOnPlane& state,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false) const = 0;

  // This interface to extrapolateToLine is intended to resemble the
  // interface of GFAbsTrackRep in old versions of genfit and is
  // implemented by default via the preceding function.
  virtual double extrapolateToLine(StateOnPlane& state,
      const TVector3& point1,
      const TVector3& point2,
      TVector3& poca,
      TVector3& dirInPoca,
      TVector3& poca_onwire,
      bool stopAtBoundary = false) const {
    TVector3 wireDir(point2 - point1);
    wireDir.Unit();
    double retval = this->extrapolateToLine(state, point1, wireDir, stopAtBoundary);
    poca = this->getPos(state);
    dirInPoca = this->getMom(state);
    dirInPoca.Unit();

    poca_onwire = point1 + wireDir*((poca - point1)*wireDir);
    
    return retval;
  }

  virtual double extrapolateToPoint(StateOnPlane& state,
      const TVector3& point,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToCylinder(StateOnPlane& state,
      double radius,
      const TVector3& linePoint = TVector3(0.,0.,0.),
      const TVector3& lineDirection = TVector3(0.,0.,1.),
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToSphere(StateOnPlane& state,
      double radius,
      const TVector3& point = TVector3(0.,0.,0.),
      bool stopAtBoundary = false) const = 0;

  /**
   * Use the Material information stored in the #TrackPoints
   */
  //virtual double extrapolateToTrackPoint() const;


  virtual unsigned int getDim() const = 0;

  virtual TVector3 getPos(const StateOnPlane& state) const = 0;

  virtual TVector3 getMom(const StateOnPlane& state) const = 0;
  TVector3 getDir(const StateOnPlane& state) const {return getMom(state).Unit();}
  virtual void getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const = 0;
  void getPosDir(const StateOnPlane& state, TVector3& pos, TVector3& dir) const {getPosMom(state, pos, dir); dir.SetMag(1.);}
  virtual TVectorD get6DState(const StateOnPlane& state) const;

  /** Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance */
  virtual void getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const = 0;
  virtual void get6DStateCov(const MeasuredStateOnPlane& state, TVectorD& stateVec, TMatrixDSym& cov) const;

  //! get the magnitude of the momentum in GeV
  virtual double getMomMag(const StateOnPlane& state) = 0;
  /** get the variance of the absolute value of the momentum  */
  virtual double getMomVar(const MeasuredStateOnPlane& state) = 0;

  int getPDG() const {return pdgCode_;}
  virtual double getCharge(const StateOnPlane& state) const = 0;
  //! get charge over momentum
  virtual double getQop(const StateOnPlane& state) const = 0;
  double getMass(const StateOnPlane& state) const;
  char getPropDir() const {return propDir_;}

  /** Get the jacobian and noise matrix of the last extrapolation  */
  virtual void getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const = 0;

  /** Get the jacobian and noise matrix of the last extrapolation if it would have been done in opposite direction  */
  virtual void getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const = 0;

  /** Get the radiation length of the material crossed in the last extrapolation.  */
  virtual double getRadiationLenght() const = 0;

  /** Calculate Jacobian of transportation numerically. Slow but accurate. Can be used to validate (semi)analytic calculations. */
  void calcJacobianNumerically(const genfit::StateOnPlane& origState,
                                   const genfit::SharedPlanePtr destPlane,
                                   TMatrixD& jacobian);

  virtual void setPosMom(StateOnPlane& state, const TVector3& pos, const TVector3& mom) const = 0;
  virtual void setPosMom(StateOnPlane& state, const TVectorD& state6) const = 0;
  virtual void setPosMomErr(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const = 0;
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const = 0;
  virtual void setPosMomCov(MeasuredStateOnPlane& state, const TVectorD& state6, const TMatrixDSym& cov6x6) const = 0;

  //! Set the sign of the charge according to charge
  virtual void setChargeSign(StateOnPlane& state, double charge) const = 0;
  virtual void setQop(StateOnPlane& state, double qop) const = 0;

  //! Set propagation direction. (-1, 0, 1) -> (backward, auto, forward)
  void setPropDir(int dir) {
    if (dir>0) propDir_ = 1;
    else if (dir<0) propDir_ = -1;
    else propDir_ = 0;
  };

  //! Switch propagation direction. Has no effect if propDir_ is set to 0.
  void switchPropDir(){propDir_ = -1*propDir_;}

  //! check if other is of same type (e.g. RKTrackRep).
  virtual bool isSameType(const AbsTrackRep* other) = 0;

  //! check if other is of same type (e.g. RKTrackRep) and has same pdg code.
  virtual bool isSame(const AbsTrackRep* other) = 0;

  virtual void Print(const Option_t* = "") const;

 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsTrackRep(const AbsTrackRep&); // copy constructor
  AbsTrackRep& operator=(const AbsTrackRep&); // assignment operator


  int pdgCode_;
  char propDir_;  // propagation direction (-1, 0, 1) -> (backward, auto, forward)

  ClassDef(AbsTrackRep,1);
};

} /* End of namespace genfit */

#endif // genfit_AbsTrackRep_h
