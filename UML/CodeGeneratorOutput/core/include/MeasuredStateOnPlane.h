#ifndef genfit_MeasuredStateOnPlane_h
#define genfit_MeasuredStateOnPlane_h

#include <TMatrixDSym.h>

#include "StateOnPlane.h"


namespace genfit {


  /** 
   *  Additional covariance matrix.
   */
class MeasuredStateOnPlane : public StateOnPlane {

 public:

  MeasuredStateOnPlane();
  MeasuredStateOnPlane(const TVectorD& state, const TMatrixDSym& cov, DetPlane* plane, AbsTrackRep* rep);
  MeasuredStateOnPlane(const StateOnPlane& state, const TMatrixDSym& cov);

  const TMatrixDSym& getCov() const {return cov_;}

  void setStateCov(const TVectorD& state, const TMatrixDSym& cov) {state_ = state; cov_ = cov;}
  void setStateCovPlane(const TVectorD& state, const TMatrixDSym& cov, DetPlane* plane) {setStatePlane(state, plane); cov_ = cov;}

 protected:

  TMatrixDSym cov_;

};

} /* End of namespace genfit */

#endif // genfit_MeasuredStateOnPlane_h
