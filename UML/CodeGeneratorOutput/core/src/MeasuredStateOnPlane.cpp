#include "MeasuredStateOnPlane.h"

namespace genfit {


  MeasuredStateOnPlane::MeasuredStateOnPlane() :
    StateOnPlane(), cov_(0,0)
  {
    ;
  }

  MeasuredStateOnPlane::MeasuredStateOnPlane(const TVectorD& state, const TMatrixDSym& cov, DetPlane* plane, AbsTrackRep* rep) :
    StateOnPlane(state, plane, rep), cov_(cov)
  {
    ;
  }

  MeasuredStateOnPlane::MeasuredStateOnPlane(const StateOnPlane& state, const TMatrixDSym& cov) :
    StateOnPlane(state), cov_(cov)
  {
    ;
  }

} /* End of namespace genfit */
