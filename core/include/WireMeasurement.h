#ifndef genfit_WireMeasurement_h
#define genfit_WireMeasurement_h

#include "AbsMeasurement.h"
#include "MeasurementOnPlane.h"


namespace genfit {

/** @brief Class for measurements in wire detectors (Straw tubes and drift chambers)
 *  which do not measure the coordinate along the wire.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Lia Lavezzi (INFN Pavia, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This hit class is not valid for any kind of plane orientation
 * choice: to use it you MUST choose a plane described by u
 * and v axes with v coincident with the wire (and u orthogonal
 * to it, obviously).
 * The hit will be described by 7 coordinates:
 * w_x1, w_y1, w_z1, w_x2, w_y2, w_z2, rdrift
 * where w_ji (with j = x, y, z and i = 1, 2) are the wire
 * extremities coordinates; rdrift = distance from the wire (u
 * coordinate in the plane)
 *
 */
class WireMeasurement : public AbsMeasurement {

 public:
  WireMeasurement(int nDim = 7);
  WireMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual SharedPlanePtr constructPlane(const StateOnPlane* state) const override;

  virtual MeasurementOnPlane constructMeasurementOnPlane(const AbsTrackRep*, const SharedPlanePtr) const override;

  virtual const TMatrixD& getHMatrix(const AbsTrackRep*) const override;

  void setMaxDistance(double d){maxDistance_ = d;}
  /**
   * select how to resolve the left/right ambiguity:
   * -1: negative (left) side on vector (track direction) x (wire direction)
   * 0: auto select (take side with smallest distance to track)
   * 1: positive (right) side on vector (track direction) x (wire direction)
   */
  void setLeftRightResolution(int lr);

  double getMaxDistance(){return maxDistance_;}
  int getLeftRightResolution() const {return leftRight_;}

 protected:

  double maxDistance_;
  double leftRight_;
};

} /* End of namespace genfit */

#endif // genfit_WireMeasurement_h
