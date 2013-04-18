#ifndef genfit_PlanarMeasurement_h
#define genfit_PlanarMeasurement_h

#include "AbsMeasurement.h"


namespace genfit {

class PlanarMeasurement : public AbsMeasurement {


 protected:
  int planeId_;
};

} /* End of namespace genfit */

#endif // genfit_PlanarMeasurement_h
