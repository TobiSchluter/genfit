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

#include <cassert>
#include <iostream>

#include "Exception.h"
#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"

namespace genfit {

KalmanFitterInfo::KalmanFitterInfo() :
  AbsFitterInfo()
{
  ;
}

KalmanFitterInfo::KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)  :
  AbsFitterInfo(trackPoint, rep)
{
  ;
}

KalmanFitterInfo::~KalmanFitterInfo() {
  ;
}


KalmanFitterInfo* KalmanFitterInfo::clone() const {
  KalmanFitterInfo* retVal = new KalmanFitterInfo(this->getTrackPoint(), this->getRep());
  if (hasReferenceState())
    retVal->setReferenceState(new ReferenceStateOnPlane(*getReferenceState()));
  if (hasForwardPrediction())
    retVal->setForwardPrediction(new MeasuredStateOnPlane(*getForwardPrediction()));
  if (hasForwardUpdate())
    retVal->setForwardUpdate(new KalmanFittedStateOnPlane(*getForwardUpdate()));
  if (hasBackwardPrediction())
    retVal->setBackwardPrediction(new MeasuredStateOnPlane(*getBackwardPrediction()));
  if (hasBackwardUpdate())
    retVal->setBackwardUpdate(new KalmanFittedStateOnPlane(*getBackwardUpdate()));

  retVal->measurementsOnPlane_.reserve(getNumMeasurements());
  for (std::vector<MeasurementOnPlane*>::const_iterator it = this->measurementsOnPlane_.begin(); it != this->measurementsOnPlane_.end(); ++it) {
    retVal->addMeasurementOnPlane(new MeasurementOnPlane(**it));
  }

  return retVal;
}


std::vector< genfit::MeasurementOnPlane* > KalmanFitterInfo::getMeasurementsOnPlane() const {
  std::vector< genfit::MeasurementOnPlane* > retVal;
  retVal.reserve(measurementsOnPlane_.size());

  for (std::vector<MeasurementOnPlane*>::const_iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    retVal.push_back(*it);
  }

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

    return calcAverageState(forwardUpdate_.get(), backwardPrediction_.get());
  }
  else { // unbiased
    if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(-1) == this->getTrackPoint()) // last measurement
      return MeasuredStateOnPlane(*forwardPrediction_);
    else if (this->getTrackPoint()->getTrack()->getPointWithMeasurement(0) == this->getTrackPoint()) // first measurement
      return MeasuredStateOnPlane(*backwardPrediction_);

    return calcAverageState(forwardPrediction_.get(), backwardPrediction_.get());
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


void KalmanFitterInfo::setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane) {
  measurementsOnPlane_.clear();

  for (std::vector<MeasurementOnPlane*>::const_iterator m = measurementsOnPlane.begin(), mend = measurementsOnPlane.end(); m < mend; ++m) {
    addMeasurementOnPlane(*m);
  }
}


void KalmanFitterInfo::setRep(const AbsTrackRep* rep) {
  rep_ = rep;

  if (referenceState_)
    referenceState_->setRep(rep);

  if (forwardPrediction_)
    forwardPrediction_->setRep(rep);

  if (forwardUpdate_)
    forwardUpdate_->setRep(rep);

  if (backwardPrediction_)
    backwardPrediction_->setRep(rep);

  if (backwardUpdate_)
    backwardUpdate_->setRep(rep);

  for (std::vector<MeasurementOnPlane*>::iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    (*it)->setRep(rep);
  }
}


void KalmanFitterInfo::deleteForwardInfo() {
  forwardPrediction_.reset();
  forwardUpdate_.reset();
}

void KalmanFitterInfo::deleteBackwardInfo() {
  backwardPrediction_.reset();
  backwardUpdate_.reset();
}

void KalmanFitterInfo::deleteReferenceInfo() {
  referenceState_.reset();
}

void KalmanFitterInfo::deleteMeasurementInfo() {
  measurementsOnPlane_.clear();
}


MeasuredStateOnPlane KalmanFitterInfo::calcAverageState(const MeasuredStateOnPlane* forwardState, const MeasuredStateOnPlane* backwardState) const {
  if (forwardState == NULL || backwardState == NULL) {
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

  if (referenceState_) {
    std::cout << "Reference state: "; referenceState_->Print();
  }
  if (forwardPrediction_) {
    std::cout << "Forward prediction_: "; forwardPrediction_->Print();
  }
  if (forwardUpdate_) {
    std::cout << "Forward update: "; forwardUpdate_->Print();
  }
  if (backwardPrediction_) {
    std::cout << "Backward prediction_: "; backwardPrediction_->Print();
  }
  if (backwardUpdate_) {
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
  /*if (!referenceState_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is NULL" << std::endl;
    return false;
  }*/

  SharedPlanePtr plane;
  if (referenceState_) {
    plane = referenceState_->getPlane();
  }
  else if (measurementsOnPlane_.size() > 0) {
    plane = measurementsOnPlane_[0]->getPlane();
  }
  else if (forwardUpdate_ || backwardUpdate_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): update w/o prediction or measurement" << std::endl;
    return false;
  }

  // see if everything else is defined wrt this plane and rep_
  if (referenceState_) {
    if (referenceState_->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): referenceState_ is not defined with the correct TrackRep" << std::endl;
      return false;
    }
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

  for (std::vector<MeasurementOnPlane*>::const_iterator it = measurementsOnPlane_.begin(); it != measurementsOnPlane_.end(); ++it) {
    if((*it)->getPlane() != plane) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct plane" << std::endl;
      return false;
    }
    if((*it)->getRep() != rep_) {
      std::cerr << "KalmanFitterInfo::checkConsistency(): measurement is not defined with the correct TrackRep" << std::endl;
      return false;
    }
  }

  // see if there is an update w/o prediction or measurement
  if (forwardUpdate_ && !forwardPrediction_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o forwardPrediction_" << std::endl;
    return false;
  }

  if (forwardUpdate_ && measurementsOnPlane_.size() == 0) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): forwardUpdate_ w/o measurement" << std::endl;
    return false;
  }


  if (backwardUpdate_ && !backwardPrediction_) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o backwardPrediction_" << std::endl;
    return false;
  }

  if (backwardUpdate_ && measurementsOnPlane_.size() == 0) {
    std::cerr << "KalmanFitterInfo::checkConsistency(): backwardUpdate_ w/o measurement" << std::endl;
    return false;
  }


  return true;
}


} /* End of namespace genfit */
