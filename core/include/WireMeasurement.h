#ifndef genfit_WireMeasurement_h
#define genfit_WireMeasurement_h

#include "AbsMeasurement.h"


namespace genfit {

class WireMeasurement : public AbsMeasurement {

 public:
  WireMeasurement() {}
  WireMeasurement(int nDim) : AbsMeasurement(nDim) {}

 protected:
  double leftRight_;
};

} /* End of namespace genfit */

#endif // genfit_WireMeasurement_h
