/* Copyright 2008-2010, Technische Universitaet Muenchen,
             2014-2015, Ludwig-Maximimilians-Universität München
   Authors: Tobias Schlüter

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_WireTimeMeasurement_h
#define genfit_WireTimeMeasurement_h

#include "AbsMeasurement.h"
#include "AbsHMatrix.h"
#include "MeasurementOnPlane.h"
#include "IDrawableMeasurement.h"


namespace genfit {

/** @brief Class for measurements in wire detectors (Straw tubes and drift chambers)
 *  which do not measure the coordinate along the wire.
 *
 *  @author Tobias Schlüter
 *
 * This is similar to WireMeasurement, but since WireMeasurement
 * stores a 7x7 covariance matrix for what is a one-dimensional
 * measurement, this class is preferable.  Protected inheritance of
 * rawHitCoords_ and rawHitCov_ makes it impossible to rewrite
 * WireMeasurement, as subclasses will access these members.
 *
 * This hit class is not valid for arbitrary choices of plane
 * orientation: to use it you MUST choose a plane described by u
 * and v axes with v coincident with the wire (and u orthogonal
 * to it, obviously).
 * The hit will be described by 7 coordinates:
 * w_x1, w_y1, w_z1, w_x2, w_y2, w_z2, rdrift
 * where w_ji (with j = x, y, z and i = 1, 2) are the wire
 * extremities coordinates; rdrift = distance from the wire (u
 * coordinate in the plane)
 *
 */
class WireTimeMeasurement : public AbsMeasurement, public display::IDrawableMeasurement {

 public:
  WireTimeMeasurement();
  WireTimeMeasurement(double TMeasured, double sigmaT, const TVector3& endPoint1, const TVector3& endPoint2, double vDrift, double vSignal, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~WireTimeMeasurement() {;}

  virtual WireTimeMeasurement* clone() const {return new WireTimeMeasurement(*this);}

  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const;

  /**  Hits with a small drift distance get a higher weight, whereas hits with
    * big drift distances become weighted down.
    * When these initial weights are used by the DAF, the smoothed track will be closer to the real
    * trajectory than if both sides are weighted with 0.5 regardless of the drift distance.
    * This helps a lot when resolving l/r ambiguities with the DAF.
    * The idea is that for the first iteration of the DAF, the wire positions are taken.
    * For small drift radii, the wire position does not bend the fit away from the
    * trajectory, whereas the wire position for hits with large drift radii is further away
    * from the trajectory and will therefore bias the fit if not weighted down.
    */
  virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const;

  virtual const AbsHMatrix* constructHMatrix(const AbsTrackRep*) const { return new HMatrix(vDrift_, vSignal_); }

  /** Reset the wire end points.
   */
  void setWireEndPoints(const TVector3& endPoint1, const TVector3& endPoint2);

  /** Set maximum drift distance. This is used to calculate the start weights of the two
   * measurementsOnPlane.
   */
  void setMaxDistance(double d){maxDistance_ = d;}
  /**
   * select how to resolve the left/right ambiguity:
   * -1: negative (left) side on vector (wire direction) x (track direction)
   * 0: mirrors enter with same weight, DAF will decide.
   * 1: positive (right) side on vector (wire direction) x (track direction)
   * where the wire direction is pointing from endPoint1 to endPoint2
   */
  void setLeftRightResolution(int lr);

  virtual bool isLeftRigthMeasurement() const {return true;}
  double getMaxDistance(){return maxDistance_;}
  int getLeftRightResolution() const {return leftRight_;}

  void drawMeasurement(TEveElementList* list, const MeasuredStateOnPlane& fittedState) const;
  void drawDetector(TEveElementList* list) const;

 private:

  double wireEndPoint1_[3]; //! Wire end point 1 (X, Y, Z)
  double wireEndPoint2_[3]; //! Wire end point 2 (X, Y, Z)
  double maxDistance_;
  double leftRight_;

  double TMeasured_;
  double sigmaT_;

  double vDrift_;
  double vSignal_;

  class HMatrix : public AbsHMatrix {
  public:
    HMatrix(double vDrift, double vSignal, bool left = false);
    virtual ~HMatrix() {;}

    virtual const TMatrixD& getMatrix() const { return H_; }
    virtual TVectorD Hv(const TVectorD& v) const {return getMatrix()*v;}
    virtual TMatrixD MHt(const TMatrixDSym& M) const {return TMatrixD(M, TMatrixD::kMultTranspose, getMatrix());}
  virtual TMatrixD MHt(const TMatrixD& M) const {return TMatrixD(M, TMatrixD::kMultTranspose, getMatrix());}
    virtual void HMHt(TMatrixDSym& M) const {M.Similarity(getMatrix());}

    virtual AbsHMatrix* clone() const { return new HMatrix(*this); }

    bool operator==(const AbsHMatrix& other) const {return this->isEqual(other);}
    bool operator!=(const AbsHMatrix& other) const {return !(this->isEqual(other));}
    virtual bool isEqual(const AbsHMatrix& other) const;

  private:
    TMatrixD H_;
  };

 public:

  ClassDef(WireTimeMeasurement, 1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_WireTimeMeasurement_h
