#ifndef genfit_KalmanFittedStateOnPlane_h
#define genfit_KalmanFittedStateOnPlane_h

#include "MeasuredStateOnPlane.h"


namespace genfit {


  /** 
   *  Additional info produced by a Kalman filter or DAF.
   */
class KalmanFittedStateOnPlane : public MeasuredStateOnPlane {

 public:
  double chiSquareIncrement_;
  
  /** 
   *  Degrees of freedom. Needs to be a double because of DAF.
   */
  double ndf_;

};

} /* End of namespace genfit */

#endif // genfit_KalmanFittedStateOnPlane_h
