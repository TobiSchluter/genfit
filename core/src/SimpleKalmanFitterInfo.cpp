/* Copyright 2013, Ludwig-Maximilians-Universität München, Technische Universität München
   Authors: Johannes Rauch & Tobias Schlüter

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

#include "SimpleKalmanFitterInfo.h"
#include "TrackPoint.h"
#include "AbsTrackRep.h"

using namespace genfit;

SimpleKalmanFitterInfo::SimpleKalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep)
  : AbsFitterInfo(trackPoint, rep)
{
  // These don't work in initializer list???
  fwPrediction_ = nullptr;
  bwPrediction_ = nullptr;
}


SimpleKalmanFitterInfo* SimpleKalmanFitterInfo::clone() const
{
  SimpleKalmanFitterInfo* res = new SimpleKalmanFitterInfo(getTrackPoint(), getRep());
  if (fwPrediction_)
    res->fwPrediction_ = std::unique_ptr<MeasuredStateOnPlane>(new MeasuredStateOnPlane(*fwPrediction_));
  if (bwPrediction_)
    res->bwPrediction_ = std::unique_ptr<MeasuredStateOnPlane>(new MeasuredStateOnPlane(*bwPrediction_));

  return res;
}


void SimpleKalmanFitterInfo::deleteForwardInfo()
{
  // unique_ptr takes care of deletion
  fwPrediction_ = nullptr;
}

void SimpleKalmanFitterInfo::deleteBackwardInfo()
{
  // unique_ptr takes care of deletion
  bwPrediction_ = nullptr;
}

MeasurementOnPlane SimpleKalmanFitterInfo::getResidual(bool biased, unsigned int iMeasurement) const
{
  // With several measurements: what is the correct definition of the
  // unbiased residual?  Including or not including the other
  // measurements?  Since the bias is introduced by the geometry, it
  // should probably be all or nothing, but then it's not clear why
  // one differentiates between the measurements at all?

  // How's the logic supposed to run?  Obviously, either this needs to
  // do some of the Kalman algebra or we need to store the residuals
  // calculated during the Kalman algebra, but then we need to
  // calculate covariances for them.  Also, these would be forward and
  // backwards residuals, thus not suited to the idea that this
  // function is meaningful independent of the track fit algorithm in use.

  assert(!biased);  //FIXME
  assert(measurements_.size() == 1); //FIXME
  assert(fwPrediction_ && bwPrediction_);

  assert(0);  // FIXME
  return MeasurementOnPlane();
}

