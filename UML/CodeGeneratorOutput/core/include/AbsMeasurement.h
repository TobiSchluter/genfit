#ifndef genfit_AbsMeasurement_h
#define genfit_AbsMeasurement_h

#include "MeasurementOnPlane.h"


namespace genfit {
class TrackPoint;
} /* End of namespace genfit */

namespace genfit {


  /** 
   *  Contains the measurement and covariance in detector coordinates. Detector and hit ids can be used to point back to the original detector hits (clusters etc.).
   */
class AbsMeasurement {

 public:

  virtual MeasurementOnPlane constructMeasurementOnPlane();


 protected:
  TVectorD rawHitCoords_;
  TMatrixDSym rawHitCov_;
  int detId_;
  int hitId_;

 public:

  /** 
   *  Can be more than one, e.g. multiple measurements in the same Si detector, left and right measurements of a wire detector etc.
   * @element-type TrackPoint
   */
  TrackPoint* rawMeasurements_;
};

} /* End of namespace genfit */

#endif // genfit_AbsMeasurement_h
