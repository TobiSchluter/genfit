#ifndef genfit_PlanarMeasurement_h
#define genfit_PlanarMeasurement_h

#include "AbsMeasurement.h"
#include "MeasurementOnPlane.h"


namespace genfit {

class AbsTrackRep;

class PlanarMeasurement : public AbsMeasurement {

 public:
  PlanarMeasurement() {}
  PlanarMeasurement(int nDim) : AbsMeasurement(nDim) {}

 protected:
  int planeId_;
};

} /* End of namespace genfit */

#endif // genfit_PlanarMeasurement_h
