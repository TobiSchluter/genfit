#ifndef  _MeasurementOnPlane_h
#define  _MeasurementOnPlane_h

#include <TMatrixD.h>

#include "MeasuredStateOnPlane.h"



namespace genfit {

class MeasurementOnPlane : public MeasuredStateOnPlane {

 public:
  TMatrixD hMatrix_;
  double weight_;

};

} /* End of namespace   */

#endif //  _MeasurementOnPlane_h
