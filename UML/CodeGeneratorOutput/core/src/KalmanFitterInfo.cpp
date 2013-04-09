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

#include "Exception.h"
#include "KalmanFitterInfo.h"
#include "Tools.h"

namespace genfit {

KalmanFitterInfo::KalmanFitterInfo() :
  AbsFitterInfo(),
  referenceState_(nullptr),
  forwardPrediction_(nullptr),
  forwardUpdate_(nullptr),
  backwardPrediction_(nullptr),
  backwardUpdate_(nullptr)
{
  ;
}

KalmanFitterInfo::KalmanFitterInfo(AbsTrackRep* rep)  :
  AbsFitterInfo(rep),
  referenceState_(nullptr),
  forwardPrediction_(nullptr),
  forwardUpdate_(nullptr),
  backwardPrediction_(nullptr),
  backwardUpdate_(nullptr)
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


KalmanFitterInfo* KalmanFitterInfo::clone() const {
  KalmanFitterInfo* retVal = new KalmanFitterInfo(this->getRep());
  if (this->referenceState_ != nullptr)
    retVal->referenceState_ = new ReferenceStateOnPlane(*(this->referenceState_));
  if (this->forwardPrediction_ != nullptr)
    retVal->forwardPrediction_ = new MeasuredStateOnPlane(*(this->forwardPrediction_));
  if (this->forwardUpdate_ != nullptr)
    retVal->forwardUpdate_ = new KalmanFittedStateOnPlane(*(this->forwardUpdate_));
  if (this->backwardPrediction_ != nullptr)
    retVal->backwardPrediction_ = new MeasuredStateOnPlane(*(this->backwardPrediction_));
  if (this->backwardUpdate_ != nullptr)
    retVal->backwardUpdate_ = new KalmanFittedStateOnPlane(*(this->backwardUpdate_));

  return retVal;
}


MeasuredStateOnPlane KalmanFitterInfo::getSmoothedState(bool biased) const {
  // TODO: Test

  if (biased) {
    if (forwardUpdate_ != nullptr && backwardPrediction_ == nullptr && backwardUpdate_ == nullptr) // last measurement
      return MeasuredStateOnPlane(*forwardUpdate_);
    else if (backwardUpdate_ != nullptr && forwardPrediction_ == nullptr && forwardUpdate_ == nullptr) // first measurement
      return MeasuredStateOnPlane(*backwardUpdate_);

    return calcSmoothedState(forwardUpdate_, backwardPrediction_);
  }
  else { // unbiased
    if (forwardPrediction_ != nullptr && backwardPrediction_ == nullptr && backwardUpdate_ == nullptr) // last measurement
      return MeasuredStateOnPlane(*forwardPrediction_);
    else if (backwardPrediction_ != nullptr && forwardPrediction_ == nullptr && forwardUpdate_ == nullptr) // first measurement
      return MeasuredStateOnPlane(*backwardPrediction_);

    return calcSmoothedState(forwardPrediction_, backwardPrediction_);
  }

}


MeasurementOnPlane KalmanFitterInfo::getResidual(bool biased, unsigned int iMeasurement) const {
  // TODO: Test

  MeasuredStateOnPlane smoothedState = getSmoothedState(biased);
  const MeasurementOnPlane& measurement = measurementsOnPlane_.at(iMeasurement);
  sharedPlanePtr plane = measurement.getPlane();

  // check equality of planes and reps
  if(*(smoothedState.getPlane()) != *plane) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }
  if(smoothedState.getRep() != measurement.getRep()) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined wrt the same TrackRep.", __LINE__,__FILE__);
    throw e;
  }

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


void KalmanFitterInfo::deleteForwardInfo() {
  if (forwardPrediction_ != nullptr) {
    delete forwardPrediction_;
    forwardPrediction_ = nullptr;
  }
  if (forwardUpdate_ != nullptr) {
    delete forwardUpdate_;
    forwardUpdate_ = nullptr;
  }
}

void KalmanFitterInfo::deleteBackwardInfo() {
  if (backwardPrediction_ != nullptr) {
    delete backwardPrediction_;
    backwardPrediction_ = nullptr;
  }
  if (backwardUpdate_ != nullptr) {
    delete backwardUpdate_;
    backwardUpdate_ = nullptr;
  }
}

void KalmanFitterInfo::deleteReferenceInfo() {
  if (referenceState_ != nullptr) {
    delete referenceState_;
    referenceState_ = nullptr;
  }
}

void KalmanFitterInfo::deleteMeasurementInfo() {
  measurementsOnPlane_.clear();
}


MeasuredStateOnPlane KalmanFitterInfo::calcSmoothedState(const MeasuredStateOnPlane* forwardState, const MeasuredStateOnPlane* backwardState) const {
  if (forwardState == nullptr || backwardState == nullptr) {
    Exception e("KalmanFitterInfo::calcSmoothedState: forwardState or backwardState is NULL.", __LINE__,__FILE__);
    throw e;
  }
  // check if both states are defined in the same plane
  if (forwardState->getPlane() != backwardState->getPlane()) {
    Exception e("KalmanFitterInfo::calcSmoothedState: forwardState and backwardState are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }

  TMatrixDSym fCovInv, bCovInv, smoothed_cov;

  const TMatrixDSym& fCov = forwardState->getCov();
  tools::invertMatrix(fCov, fCovInv);

  const TMatrixDSym& bCov = backwardState->getCov();
  tools::invertMatrix(bCov, bCovInv);

  tools::invertMatrix(fCovInv + bCovInv, smoothed_cov);

  return MeasuredStateOnPlane(smoothed_cov * (fCovInv*forwardState->getState() + bCovInv*backwardState->getState()), // smoothed state
                              smoothed_cov,
                              forwardUpdate_->getPlane(),
                              forwardUpdate_->getRep());
}


} /* End of namespace genfit */
