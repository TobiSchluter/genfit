#ifndef genfit_WireMeasurement_h
#define genfit_WireMeasurement_h

#include "AbsMeasurement.h"
#include "MeasurementOnPlane.h"


namespace genfit {

class WireMeasurement : public AbsMeasurement {
  static const double HMatrixContent_[5];
  static const TMatrixD HMatrix_;

 public:
  WireMeasurement() {}
  WireMeasurement(int nDim) : AbsMeasurement(nDim) {}

  virtual MeasurementOnPlane constructMeasurementOnPlane(const AbsTrackRep*, const MeasuredStateOnPlane&) const override;

 protected:
  double leftRight_;
};

} /* End of namespace genfit */

#endif // genfit_WireMeasurement_h
