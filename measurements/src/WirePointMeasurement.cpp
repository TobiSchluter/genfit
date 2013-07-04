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

#include "WirePointMeasurement.h"

#include "Exception.h"
#include "RKTrackRep.h"

#include <cassert>


namespace genfit {

WirePointMeasurement::WirePointMeasurement(int nDim)
  : WireMeasurement(nDim)
{
  assert(nDim >= 8);
}

WirePointMeasurement::WirePointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : WireMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint)
{
  assert(rawHitCoords_.GetNrows() >= 8);
}


std::vector<MeasurementOnPlane*> WirePointMeasurement::constructMeasurementsOnPlane(const AbsTrackRep* rep, const SharedPlanePtr plane) const
{
  MeasurementOnPlane* mopR = new MeasurementOnPlane(TVectorD(2),
       TMatrixDSym(2),
       plane, rep, getHMatrix(rep));

  mopR->getState()(0) = rawHitCoords_(6);
  mopR->getState()(1) = rawHitCoords_(7);

  mopR->getCov()(0,0) = rawHitCov_(6,6);
  mopR->getCov()(1,0) = rawHitCov_(7,6);
  mopR->getCov()(0,1) = rawHitCov_(6,7);
  mopR->getCov()(1,1) = rawHitCov_(7,7);


  MeasurementOnPlane* mopL = new MeasurementOnPlane(*mopR);
  mopL->getState()(0) *= -1;

  // set left/right weights
  if (leftRight_ < 0) {
    mopL->setWeight(1);
    mopR->setWeight(0);
  }
  else if (leftRight_ > 0) {
    mopL->setWeight(0);
    mopR->setWeight(1);
  }

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mopL);
  retVal.push_back(mopR);
  return retVal;
}

const TMatrixD& WirePointMeasurement::getHMatrix(const AbsTrackRep* rep) const {
  if (dynamic_cast<const RKTrackRep*>(rep) != NULL) {
    static const double HMatrixContent[10] = {0, 0, 0, 1, 0,
                                              0, 0, 0, 0, 1};
    static const TMatrixT<double> HMatrix(2,5, HMatrixContent);

    return HMatrix;
  }
  else {
    Exception exc("WirePointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }
}

} /* End of namespace genfit */
