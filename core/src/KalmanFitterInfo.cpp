/* Copyright 2008-2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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

#include "KalmanFitterInfo.h"

#include <assert.h>
#include <iostream>

#include "Exception.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

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

KalmanFitterInfo::KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)  :
  AbsFitterInfo(trackPoint, rep),
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

  for (MeasurementOnPlane* m : measurementsOnPlane_) {
    if (m != nullptr)
      delete m;
  }
}


KalmanFitterInfo* KalmanFitterInfo::clone() const {
  KalmanFitterInfo* retVal = new KalmanFitterInfo(this->getTrackPoint(), this->getRep());
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

  retVal->measurementsOnPlane_.reserve(this->measurementsOnPlane_.size());
  for (MeasurementOnPlane* mop : this->measurementsOnPlane_)
    retVal->measurementsOnPlane_.push_back(new MeasurementOnPlane(*mop));

  return retVal;
}


MeasurementOnPlane KalmanFitterInfo::getAvgWeightedMeasurementOnPlane() const {

  MeasurementOnPlane retVal(*(measurementsOnPlane_[0]));
  retVal.setWeight(1);

  if(measurementsOnPlane_.size() == 1) {
    if (retVal.getWeight() != 1.)
      retVal.getCov() *= 1. / retVal.getWeight();
  }
  else { // more than one hit
    retVal.getState().Zero();
    retVal.getCov().Zero();

    TMatrixDSym covInv;
    std::vector<TMatrixDSym> weightedCovInvs;

    for(unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {

      if (i>0) {
        // make sure we have compatible measurement types
        // TODO: replace with Exceptions!
        assert(measurementsOnPlane_[i]->getPlane() == measurementsOnPlane_[0]->getPlane());
        assert(measurementsOnPlane_[i]->getHMatrix() == measurementsOnPlane_[0]->getHMatrix());
      }

      tools::invertMatrix(measurementsOnPlane_[i]->getCov(), covInv); // invert cov
      covInv *= measurementsOnPlane_[i]->getWeight(); // weigh cov
      weightedCovInvs.push_back(covInv); // cov is already inverted and weighted
      retVal.getCov() += covInv; // cov is already inverted and weighted
    }

    // invert fHitCov
    tools::invertMatrix(retVal.getCov());

    //set the weighted-mean coord
    for(unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {
      retVal.getState() += weightedCovInvs[i] * measurementsOnPlane_[i]->getState();
    }
    retVal.getState() *= retVal.getCov();
  }

  return retVal;
}


MeasuredStateOnPlane KalmanFitterInfo::getFittedState(bool biased) const {
  // TODO: Test

  if (biased) {
    if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(-1) == this->getTrackPoint()) // last measurement
      return MeasuredStateOnPlane(*forwardUpdate_);
    else if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(0) == this->getTrackPoint()) // first measurement
      return MeasuredStateOnPlane(*backwardUpdate_);

    return calcAverageState(forwardUpdate_, backwardPrediction_);
  }
  else { // unbiased
    if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(-1) == this->getTrackPoint()) // last measurement
      return MeasuredStateOnPlane(*forwardPrediction_);
    else if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(0) == this->getTrackPoint()) // first measurement
      return MeasuredStateOnPlane(*backwardPrediction_);

    return calcAverageState(forwardPrediction_, backwardPrediction_);
  }

}


MeasurementOnPlane KalmanFitterInfo::getResidual(bool biased, unsigned int iMeasurement) const {
  // TODO: Test

  MeasuredStateOnPlane smoothedState = getFittedState(biased);
  const MeasurementOnPlane* measurement = measurementsOnPlane_.at(iMeasurement);
  SharedPlanePtr plane = measurement->getPlane();

  // check equality of planes and reps
  if(*(smoothedState.getPlane()) != *plane) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined in the same plane.", __LINE__,__FILE__);
    throw e;
  }
  if(smoothedState.getRep() != measurement->getRep()) {
    Exception e("KalmanFitterInfo::getResidual: smoothedState and measurement are not defined wrt the same TrackRep.", __LINE__,__FILE__);
    throw e;
  }

  const TMatrixD& H = measurement->getHMatrix();

  TVectorD res = measurement->getState() - (H * smoothedState.getState());

  TMatrixDSym cov(smoothedState.getCov());
  cov.Similarity(H);
  cov += measurement->getCov();

  return MeasurementOnPlane(res, cov, plane, smoothedState.getRep(), H, measurement->getWeight());
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


void KalmanFitterInfo::setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane) {
  for (MeasurementOnPlane* m : measurementsOnPlane_) {
    if (m != nullptr)
      delete m;
  }

  measurementsOnPlane_ = measurementsOnPlane;
}


void KalmanFitterInfo::setRep(const AbsTrackRep* rep) {
  rep_ = rep;

  if (referenceState_ != nullptr)
    referenceState_->setRep(rep);

  if (forwardPrediction_ != nullptr)
    forwardPrediction_->setRep(rep);

  if (forwardUpdate_ != nullptr)
    forwardUpdate_->setRep(rep);

  if (backwardPrediction_ != nullptr)
    backwardPrediction_->setRep(rep);

  if (backwardUpdate_ != nullptr)
    backwardUpdate_->setRep(rep);

  for (MeasurementOnPlane* m : measurementsOnPlane_)
    m->setRep(rep);
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


MeasuredStateOnPlane KalmanFitterInfo::calcAverageState(const MeasuredStateOnPlane* forwardState, const MeasuredStateOnPlane* backwardState) const {
  if (forwardState == nullptr || backwardState == nullptr) {
    Exception e("KalmanFitterInfo::calcAverageState: forwardState or backwardState is NULL.", __LINE__,__FILE__);
    throw e;
  }
  // check if both states are defined in the same plane
  if (forwardState->getPlane() != backwardState->getPlane()) {
    Exception e("KalmanFitterInfo::calcAverageState: forwardState and backwardState are not defined in the same plane.", __LINE__,__FILE__);
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


void KalmanFitterInfo::Print(const Option_t*) const {
  std::cout << "genfit::KalmanFitterInfo. Belongs to TrackPoint " << trackPoint_ << "; TrackRep " <<  rep_ << "; statusFlag = " << statusFlag_ << "\n";

  for (unsigned int i=0; i<measurementsOnPlane_.size(); ++i) {
    std::cout << "MeasurementOnPlane Nr " << i <<": "; measurementsOnPlane_[i]->Print();
  }

  if (referenceState_ != nullptr) {
    std::cout << "Reference state: "; referenceState_->Print();
  }
  if (forwardPrediction_ != nullptr) {
    std::cout << "Forward prediction_: "; forwardPrediction_->Print();
  }
  if (forwardUpdate_ != nullptr) {
    std::cout << "Forward update: "; forwardUpdate_->Print();
  }
  if (backwardPrediction_ != nullptr) {
    std::cout << "Backward prediction_: "; backwardPrediction_->Print();
  }
  if (backwardUpdate_ != nullptr) {
    std::cout << "Backward update: "; backwardUpdate_->Print();
  }

}


bool KalmanFitterInfo::checkConsistency() const {

  // check if in a TrackPoint
  if (!trackPoint_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): trackPoint_ is NULL" << std::endl;
    return false;
  }

  // check if there is a reference state
  if (!referenceState_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is NULL" << std::endl;
    return false;
  }

  SharedPlanePtr plane = referenceState_->getPlane();

  // see if everything else is defined wrt this plane and rep_
  if (referenceState_->getRep() != rep_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is not defined with the correct TrackRep" << std::endl;
    return false;
  }

  if (forwardPrediction_) {
    if(forwardPrediction_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(forwardPrediction_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardPrediction_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
  }
  if (forwardUpdate_) {
    if(forwardUpdate_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(forwardUpdate_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
  }

  if (backwardPrediction_) {
    if(backwardPrediction_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(backwardPrediction_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardPrediction_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
  }
  if (backwardUpdate_) {
    if(backwardUpdate_->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ is not defined with the correct plane" << std::endl;
      return false;
    }
    if(backwardUpdate_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
  }

  for (MeasurementOnPlane* m : measurementsOnPlane_){
    if(m->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct plane" << std::endl;
      return false;
    }
    if(m->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct TrackRep" << std::endl;
      return false;
    }
  }


  // see if there is an update w/o prediction
  if (forwardUpdate_ && !forwardPrediction_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o forwardPrediction_" << std::endl;
    return false;
  }

  if (backwardUpdate_ && !backwardPrediction_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o backwardPrediction_" << std::endl;
    return false;
  }


  return true;
}


} /* End of namespace genfit */
