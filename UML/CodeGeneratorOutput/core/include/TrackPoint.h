#ifndef genfit_TrackPoint_h
#define genfit_TrackPoint_h

#include <vector>

#include "AbsMeasurement.h"
#include "KalmanFitterInfo.h"
#include "MaterialInfo.h"

namespace genfit {
class Track;
} /* End of namespace genfit */

namespace genfit {

class TrackPoint {

 public:
  double sortingParameter_;

 public:

  /**
   * @element-type Track
   */
  Track *trackPoints_;

  /** 
   *  Can be more than one, e.g. multiple measurements in the same Si detector, left and right measurements of a wire detector etc.
   * @element-type AbsMeasurement
   */
  std::vector< AbsMeasurement* > rawMeasurements_;

  /**
   * @element-type KalmanFitterInfo
   */
  std::vector< KalmanFitterInfo > fitterInfos_;

  MaterialInfo * material_;
};

} /* End of namespace genfit */

#endif // genfit_TrackPoint_h
