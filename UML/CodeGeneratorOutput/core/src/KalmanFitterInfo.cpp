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

#include "KalmanFitterInfo.h"

namespace genfit {

KalmanFitterInfo::KalmanFitterInfo() :
  referenceState_(NULL),
  forwardPrediction_(NULL),
  forwardUpdate_(NULL),
  backwardPrediction_(NULL),
  backwardUpdate_(NULL),
  rep_(NULL)
{
  ;
}

KalmanFitterInfo::KalmanFitterInfo(AbsTrackRep* rep)  :
  referenceState_(NULL),
  forwardPrediction_(NULL),
  forwardUpdate_(NULL),
  backwardPrediction_(NULL),
  backwardUpdate_(NULL),
  rep_(rep)
{
  ;
}

KalmanFitterInfo::~KalmanFitterInfo() {
  if (referenceState_ != NULL)
    delete referenceState_;
  if (forwardPrediction_ != NULL)
    delete forwardPrediction_;
  if (forwardUpdate_ != NULL)
    delete forwardUpdate_;
  if (backwardPrediction_ != NULL)
    delete backwardPrediction_;
  if (backwardUpdate_ != NULL)
    delete backwardUpdate_;
}


MeasuredStateOnPlane KalmanFitterInfo::getBiasedSmoothedState() const {
  // TODO: implement
}

MeasuredStateOnPlane KalmanFitterInfo::getUnbiasedSmoothedState() const {
  // TODO: implement
}

StateOnPlane KalmanFitterInfo::getBiasedResidual() const {
  // TODO: implement
}

StateOnPlane KalmanFitterInfo::getUnbiasedResidual() const {
  // TODO: implement
}


void KalmanFitterInfo::setReferenceState(ReferenceStateOnPlane* referenceState) {
  if (referenceState_ != NULL)
    delete referenceState_;
  referenceState_ = referenceState;
}

void KalmanFitterInfo::setForwardPrediction(MeasuredStateOnPlane* forwardPrediction) {
  if (forwardPrediction_ != NULL)
    delete forwardPrediction_;
  forwardPrediction_ = forwardPrediction;
}

void KalmanFitterInfo::setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate) {
  if (forwardUpdate_ != NULL)
    delete forwardUpdate_;
  forwardUpdate_ = forwardUpdate;
}

void KalmanFitterInfo::setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction) {
  if (backwardPrediction_ != NULL)
    delete backwardPrediction_;
  backwardPrediction_ = backwardPrediction;
}

void KalmanFitterInfo::setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate) {
  if (backwardUpdate_ != NULL)
    delete backwardUpdate_;
  backwardUpdate_ = backwardUpdate;
}


} /* End of namespace genfit */
