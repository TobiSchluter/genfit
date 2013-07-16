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

#include "PlanarMeasurement.h"

#include "Exception.h"
#include "RKTrackRep.h"

#include <cassert>


namespace genfit {

PlanarMeasurement::PlanarMeasurement(int nDim)
  : AbsMeasurement(nDim), physicalPlane_(), planeId_(-1)
{
  assert(nDim >= 1);
}

PlanarMeasurement::PlanarMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint), physicalPlane_(), planeId_(-1)
{
  assert(rawHitCoords_.GetNrows() >= 1);
}


SharedPlanePtr PlanarMeasurement::constructPlane(const StateOnPlane* state) const {
  if (!physicalPlane_) {
    Exception exc("PlanarMeasurement::constructPlane(): No plane has been set!", __LINE__,__FILE__);
    throw exc;
  }
  return physicalPlane_;
}


std::vector<MeasurementOnPlane*> PlanarMeasurement::constructMeasurementsOnPlane(const AbsTrackRep* rep, const SharedPlanePtr plane) const {

  MeasurementOnPlane* mop = new MeasurementOnPlane(rawHitCoords_,
       rawHitCov_,
       plane, rep, getHMatrix(rep));

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mop);
  return retVal;
}


const TMatrixD& PlanarMeasurement::getHMatrix(const AbsTrackRep* rep) const {

  if (dynamic_cast<const RKTrackRep*>(rep) != NULL) {
    switch(rawHitCoords_.GetNrows()) {
    case 1:
      static const double HMatrixContent1[5] = {0, 0, 0, 1, 0};
      static const TMatrixT<double> HMatrix1(1,5, HMatrixContent1);
      return HMatrix1;

    case 2:
      static const double HMatrixContent2[10] = {0, 0, 0, 1, 0,
                                                0, 0, 0, 0, 1};
      static const TMatrixT<double> HMatrix2(2,5, HMatrixContent2);
      return HMatrix2;

    default:
      Exception exc("SpacepointMeasurement default implementation can only handle 1D (strip) or 2D (pixel) measurements!", __LINE__,__FILE__);
      throw exc;
    }
  }
  else {
    Exception exc("SpacepointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }
}

void PlanarMeasurement::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::PlanarMeasurement.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::PlanarMeasurement thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      DetPlane* p = 0;
      R__b >> p;
      this->physicalPlane_.reset(p);
      R__b >> planeId_;
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      R__b << physicalPlane_.get();
      R__b << planeId_;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
