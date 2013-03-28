#ifndef genfit_MeasuredStateOnPlane_h
#define genfit_MeasuredStateOnPlane_h

#include <TMatrixDSym.h>

#include "StateOnPlane.h"


namespace genfit {


  /** 
   *  Additional covariance matrix.
   */
class MeasuredStateOnPlane : public StateOnPlane {


 protected:
  TMatrixDSym cov_;

};

} /* End of namespace genfit */

#endif // genfit_MeasuredStateOnPlane_h
