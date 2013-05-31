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
/** @addtogroup genfit
 * @{
 */

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


  /** 
   *  This class collects all information needed and produced by a Kalman filter or DAF and is specific to one #GFAbsTrackRep of the #GFTrack.
   */
class KalmanFitterInfo : public AbsFitterInfo {

 public:

  KalmanFitterInfo();
  KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);
  ~KalmanFitterInfo();

  virtual KalmanFitterInfo* clone() const override;

  const ReferenceStateOnPlane* getReferenceState() const {return referenceState_;}
  const MeasuredStateOnPlane* getForwardPrediction() const {return forwardPrediction_;}
  const KalmanFittedStateOnPlane* getForwardUpdate() const {return forwardUpdate_;}
  const MeasuredStateOnPlane* getBackwardPrediction() const {return backwardPrediction_;}
  const KalmanFittedStateOnPlane* getBackwardUpdate() const {return backwardUpdate_;}
  const std::vector< genfit::MeasurementOnPlane >& getMeasurementsOnPlane() const {return measurementsOnPlane_;}
  const MeasurementOnPlane& getMeasurementOnPlane(unsigned int i) const {return measurementsOnPlane_.at(i);}

  bool hasReferenceState() const {return (referenceState_ != nullptr);}
  bool hasForwardPrediction() const {return (forwardPrediction_ != nullptr);}
  bool hasForwardUpdate() const {return (forwardUpdate_ != nullptr);}
  bool hasBackwardPrediction() const {return (backwardPrediction_ != nullptr);}
  bool hasBackwardUpdate() const {return (backwardUpdate_ != nullptr);}
  unsigned int getNumMeasurements() {return measurementsOnPlane_.size();}

  /** Get unbiased (default) or biased smoothed state
   */
  MeasuredStateOnPlane getSmoothedState(bool biased = false) const;
  /** Get unbiased (default) or biased residual from ith measurement
   */
  MeasurementOnPlane getResidual(bool biased = false, unsigned int iMeasurement = 0) const override; // also calculates covariance of the residual

  void setReferenceState(ReferenceStateOnPlane* referenceState);
  void setForwardPrediction(MeasuredStateOnPlane* forwardPrediction);
  void setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate);
  void setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction);
  void setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate);
  void setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane >& measurementsOnPlane) {measurementsOnPlane_ = measurementsOnPlane;}

  void deleteForwardInfo() override;
  void deleteBackwardInfo() override;
  void deleteReferenceInfo() override;
  void deleteMeasurementInfo() override;

 private:

  MeasuredStateOnPlane calcAverageState(const MeasuredStateOnPlane* forwardState, const MeasuredStateOnPlane* backwardState) const;

  ReferenceStateOnPlane* referenceState_; // Ownership
  MeasuredStateOnPlane* forwardPrediction_; // Ownership
  KalmanFittedStateOnPlane* forwardUpdate_; // Ownership
  MeasuredStateOnPlane* backwardPrediction_; // Ownership
  KalmanFittedStateOnPlane* backwardUpdate_; // Ownership

  /** 
   *  Number of measurements must be equal to size of #fRawMeasurements in #GFTrackPoint.
   * @element-type MeasurementOnPlane
   */
  std::vector< genfit::MeasurementOnPlane > measurementsOnPlane_;


  //ClassDef(KalmanFitterInfo,1)

};

} /* End of namespace genfit */

#endif // genfit_KalmanFitterInfo_h
