/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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
/**
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */


/** @addtogroup utilities
 * @{
 */

#ifndef genfit_MeasurementOnPlaneCreator_h
#define genfit_MeasurementOnPlaneCreator_h

#include <AbsMeasurement.h>
#include <HelixTrackModel.h>

#include <TObject.h>
#include <TVector3.h>


namespace genfit {

/** @brief Detector plane genfit geometry class
 *
 * A detector plane is the principle object to define coordinate systems for
 * track fitting in genfit. Since a particle trajectory is a
 * one-dimensional object (regardless of any specific parameterization)
 * positions with respect to the track are always measured in a plane.
 *
 * Which plane is chosen depends on the type of detector. Fixed plane
 * detectors have their detector plane defined by their mechanical setup. While
 * wire chambers or time projection chambers might want to define a detector
 * plane more flexibly.
 *
 * This class parameterizes a plane in terms of an origin vector o
 * and two plane-spanning directions u and v.
 */

enum eMeasurementType { Pixel = 0,
        Spacepoint,
        ProlateSpacepoint,
        StripU,
        StripV,
        Wire,
        WirePoint };

class MeasurementCreator : public TObject {

  /**
   * Helix track model for testing purposes
   */

 public:

  // Constructors/Destructors ---------
  MeasurementCreator();

  ~MeasurementCreator() {delete trackModel_;}

  void setTrackModel(const HelixTrackModel* model) {delete trackModel_; trackModel_ = model;}
  void setResolution(double resolution) {resolution_ = resolution;}
  void setResolutionWire(double resolutionWire) {resolutionWire_ = resolutionWire;}
  void setOutlierProb(double outlierProb) {outlierProb_ = outlierProb;}
  void setOutlierRange(double outlierRange) {outlierRange_ = outlierRange;}
  void setThetaDetPlane(double thetaDetPlane) {thetaDetPlane_ = thetaDetPlane;}
  void setPhiDetPlane(double phiDetPlane) {phiDetPlane_ = phiDetPlane;}
  void setWireDir(const TVector3 wireDir) {wireDir_ = wireDir; wireDir_.SetMag(1.);}
  void setMinDrift(double minDrift) {minDrift_ = minDrift;}
  void setMaxDrift(double maxDrift) {maxDrift_ = maxDrift;}
  void setIdealLRResolution(bool idealLRResolution) {idealLRResolution_ = idealLRResolution;}
  void setUseSkew(bool useSkew) {useSkew_ = useSkew;}
  void setSkewAngle(double skewAngle) {skewAngle_ = skewAngle;}
  void setNSuperLayer(int nSuperLayer) {nSuperLayer_ = nSuperLayer;}
  void setDebug(bool debug) {debug_ = debug;}


  AbsMeasurement* create(eMeasurementType, double tracklength, bool& outlier, int& lr);
  AbsMeasurement* create(eMeasurementType type, double tracklength) {
    bool dummy1;
    int dummy2;
    return create(type, tracklength, dummy1, dummy2);
  }

  void reset();

 private:

  const HelixTrackModel* trackModel_; // ownership

  double resolution_;  // cm; resolution of generated measurements
  double resolutionWire_;  // cm; resolution in wire direction of generated measurements (wire and prolate sp measurements)

  double outlierProb_;
  double outlierRange_;

  // planarMeasurement specific
  double thetaDetPlane_; // degree
  double phiDetPlane_; // degree

  // WireMeasurement specific
  int wireCounter_;
  TVector3 wireDir_;
  double minDrift_;
  double maxDrift_;
  bool idealLRResolution_; // resolve the l/r ambiguities of the wire measurements
  bool useSkew_;
  double skewAngle_;
  int nSuperLayer_;

  // misc
  int measurementCounter_;
  bool debug_;


 public:
  ClassDef(MeasurementCreator,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_MeasurementOnPlaneCreator_h
