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

  virtual ~AbsTrackRep() = 0;

  virtual AbsTrackRep* clone() const = 0;

  /** Extrapolates the stateInput to plane, and returns the extrapolation length
   * and, via reference, the extrapolated statePrediction.
   * If stopAtBoundary is true, the extrapolation stops as soon as a material boundary is encountered.
   */
  virtual double extrapolateToPlane(const StateOnPlane& stateInput,
      StateOnPlane& statePrediction,
      sharedPlanePtr plane,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToLine(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToPoint(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& point,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToCylinder(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& linePoint,
      const TVector3& lineDirection,
      double radius,
      bool stopAtBoundary = false) const = 0;

  virtual double extrapolateToSphere(const StateOnPlane* stateInput,
      StateOnPlane* statePrediction,
      const TVector3& point,
      double radius,
      bool stopAtBoundary = false) const = 0;

  virtual TVector3 getPos(const StateOnPlane* stateInput) const = 0;

  virtual TVector3 getMom(const StateOnPlane* stateInput) const = 0;
  virtual void getPosMom(const StateOnPlane* stateInput, TVector3& pos, TVector3& mom) const = 0;

  /**
   * Translates MeasuredStateOnPlane into 3D position, momentum and 6x6 covariance
   */
  virtual void getPosMomCov(const MeasuredStateOnPlane* stateInput, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const = 0;


  int getPDG() const {return pdgCode_;}
  virtual double getCharge() const = 0;

  /**
   * Use the Material information stored in the #TrackPoints
   */
  //virtual double extrapolateToTrackPoint() const;


 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsTrackRep(const AbsTrackRep&) = default; // copy constructor
  AbsTrackRep& operator=(const AbsTrackRep&) = default; // assignment operator


  int pdgCode_;

};

} /* End of namespace genfit */

#endif // genfit_AbsTrackRep_h
