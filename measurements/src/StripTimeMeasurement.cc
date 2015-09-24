/* Copyright 2008-2015, Technische Universitaet Muenchen,
                        Ludwig-Maximilians-Universität München
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch
            & Tobias Schlüter

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

#include "StripTimeMeasurement.h"

#include <Exception.h>
#include <RKTrackRep.h>
#include <RKTrackRepEnergy.h>
#include <RKTrackRepTime.h>
#include <HMatrixUT6.h>

#include <TEveBox.h>

#include <cassert>


namespace genfit {

StripTimeMeasurement::StripTimeMeasurement()
  : AbsMeasurement(2), physicalPlane_(), planeId_(-1), stripV_(false)
{}

StripTimeMeasurement::StripTimeMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint), physicalPlane_(), planeId_(-1), stripV_(false)
{
  assert(rawHitCoords_.GetNrows() == 2);
}


SharedPlanePtr StripTimeMeasurement::constructPlane(const StateOnPlane&) const {
  if (!physicalPlane_) {
    Exception exc("StripTimeMeasurement::constructPlane(): No plane has been set!", __LINE__,__FILE__);
    throw exc;
  }
  return physicalPlane_;
}


std::vector<MeasurementOnPlane*> StripTimeMeasurement::constructMeasurementsOnPlane(const StateOnPlane& state) const {

  MeasurementOnPlane* mop = new MeasurementOnPlane(rawHitCoords_,
       rawHitCov_,
       state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mop);
  return retVal;
}


const AbsHMatrix* StripTimeMeasurement::constructHMatrix(const AbsTrackRep* rep) const {

  if (dynamic_cast<const RKTrackRepTime*>(rep)) {
    if (stripV_)
      ; //return new HMatrixVT6();
    else
      return new HMatrixUT6();
  }

  Exception exc("SpacepointMeasurement default implementation can only handle state vectors of type RKTrackRepTime!", __LINE__,__FILE__);
  throw exc;
}

void StripTimeMeasurement::drawDetector(TEveElementList *list) const
{
  SharedPlanePtr physPlane = this->getPhysicalPlane();
  const TVector3& o(physPlane->getO());
  const TVector3& u(physPlane->getU());
  const TVector3& v(physPlane->getV());
  TEveBox* box = genfit::display::boxCreator(o, u, v, 4, 4, 0.01);
  box->SetMainColor(kGray);
  box->SetMainTransparency(50);

  list->AddElement(box);
}

void StripTimeMeasurement::drawMeasurement(TEveElementList *list, const MeasuredStateOnPlane& fittedState) const
{
  const TVector3& track_pos = fittedState.getPos();
  const TVector3& o = fittedState.getPlane()->getO();
  const TVector3& u = fittedState.getPlane()->getU();
  const TVector3& v = fittedState.getPlane()->getV();

  double hit_u = getRawHitCoords()(0);
  const TMatrixDSym& hit_cov = getRawHitCov();

  double errorScale_ = 1;
  double plane_size = 4;

  TEveBox* hit_box;
  TVector3 stripDir3 = u;
  TVector3 stripDir3perp = -v;
  if (stripV_) {
    stripDir3 = v;
    stripDir3perp = u;
  }
  TVector3 move = stripDir3perp*(stripDir3perp*(track_pos-o));
  hit_box = genfit::display::boxCreator((o + move + hit_u*stripDir3), stripDir3, stripDir3perp, errorScale_*std::sqrt(hit_cov(0,0)), plane_size, 0.0105);
  hit_box->SetMainColor(kYellow);
  hit_box->SetMainTransparency(0);

  list->AddElement(hit_box);
}

void StripTimeMeasurement::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::StripTimeMeasurement.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::StripTimeMeasurement thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      char flag;
      R__b >> flag;
      physicalPlane_.reset();
      if (flag) {
        physicalPlane_.reset(new DetPlane());
        physicalPlane_->Streamer(R__b);
      }
      R__b >> planeId_;
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      if (physicalPlane_) {
        R__b << (char)1;
        physicalPlane_->Streamer(R__b);
      } else {
        R__b << (char)0;
      }
      R__b << planeId_;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
