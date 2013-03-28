#ifndef genfit_KalmanFitterInfo_h
#define genfit_KalmanFitterInfo_h

#include <vector>


#include "AbsFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"
#include "ReferenceStateOnPlane.h"
#include "StateOnPlane.h"

namespace genfit {
class AbsTrackRep;
} /* End of namespace genfit */

namespace genfit {


  /** 
   *  This class collects all information needed and produced by a Kalman filter or DAF and is specific to one #GFAbsTrackRep of the #GFTrack.
   */
class KalmanFitterInfo : public AbsFitterInfo {

 public:

  virtual MeasuredStateOnPlane getBiasedSmoothedState();

  virtual MeasuredStateOnPlane getUnbiasedSmoothedState();

  virtual StateOnPlane getBiasedResidual();

  virtual StateOnPlane getUnbiasedResidual();

 public:

  ReferenceStateOnPlane* referenceState_;

  MeasuredStateOnPlane* forwardPrediction_;

  KalmanFittedStateOnPlane* forwardUpdate_;

  MeasuredStateOnPlane* backwardPrediction_;

  KalmanFittedStateOnPlane* backwardUpdate_;


  /** 
   *  Number of measurements must be equal to size of #fRawMeasurements in #GFTrackPoint.
   * @element-type MeasurementOnPlane
   */
  std::vector< MeasurementOnPlane > measurementsOnPlane_;

  AbsTrackRep* rep_;
};

} /* End of namespace genfit */

#endif // genfit_KalmanFitterInfo_h
