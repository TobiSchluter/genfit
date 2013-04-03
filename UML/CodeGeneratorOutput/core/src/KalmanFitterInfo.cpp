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

#include <assert.h>

#include "KalmanFitterInfo.h"

namespace genfit {

KalmanFitterInfo::KalmanFitterInfo() :
  referenceState_(nullptr),
  forwardPrediction_(nullptr),
  forwardUpdate_(nullptr),
  backwardPrediction_(nullptr),
  backwardUpdate_(nullptr),
  rep_(nullptr)
{
  ;
}

KalmanFitterInfo::KalmanFitterInfo(AbsTrackRep* rep)  :
  referenceState_(nullptr),
  forwardPrediction_(nullptr),
  forwardUpdate_(nullptr),
  backwardPrediction_(nullptr),
  backwardUpdate_(nullptr),
  rep_(rep)
{
  ;
}

KalmanFitterInfo::~KalmanFitterInfo() {
  if (referenceState_ != nullptr)
    delete referenceState_;
  if (forwardPrediction_ != nullptr)
    delete forwardPrediction_;
  if (forwardUpdate_ != nullptr)
    delete forwardUpdate_;
  if (backwardPrediction_ != nullptr)
    delete backwardPrediction_;
  if (backwardUpdate_ != nullptr)
    delete backwardUpdate_;
}


MeasuredStateOnPlane KalmanFitterInfo::getSmoothedState(bool biased) const {
  // TODO: implement
  return MeasuredStateOnPlane();
}


MeasurementOnPlane KalmanFitterInfo::getResidual(bool biased, unsigned int iMeasurement) const {
  // TODO: Test

  MeasuredStateOnPlane smoothedState = getSmoothedState(biased);
  const MeasurementOnPlane& measurement = measurementsOnPlane_.at(iMeasurement);
  const DetPlane* plane = measurement.getPlane();

  // check equality of planes and reps
  assert(*(smoothedState.getPlane()) == *plane); // TODO: replace assertion
  assert(smoothedState.getRep() == measurement.getRep()); // TODO: replace assertion

  const TMatrixD& H = measurement.getHMatrix();

  TVectorD res = measurement.getState() - (H * smoothedState.getState());

  TMatrixDSym cov(smoothedState.getCov());
  cov.Similarity(H);
  cov += measurement.getCov();

  return MeasurementOnPlane(res, cov, plane, smoothedState.getRep(), H, measurement.getWeight());
}


void KalmanFitterInfo::setReferenceState(ReferenceStateOnPlane* referenceState) {
  if (referenceState_ != nullptr)
    delete referenceState_;
  referenceState_ = referenceState;
}

void KalmanFitterInfo::setForwardPrediction(MeasuredStateOnPlane* forwardPrediction) {
  if (forwardPrediction_ != nullptr)
    delete forwardPrediction_;
  forwardPrediction_ = forwardPrediction;
}

void KalmanFitterInfo::setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate) {
  if (forwardUpdate_ != nullptr)
    delete forwardUpdate_;
  forwardUpdate_ = forwardUpdate;
}

void KalmanFitterInfo::setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction) {
  if (backwardPrediction_ != nullptr)
    delete backwardPrediction_;
  backwardPrediction_ = backwardPrediction;
}

void KalmanFitterInfo::setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate) {
  if (backwardUpdate_ != nullptr)
    delete backwardUpdate_;
  backwardUpdate_ = backwardUpdate;
}


} /* End of namespace genfit */
