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

#include <Exception.h>
#include <RKTrackRep.h>
#include <RKTrackRepEnergy.h>
#include <RKTrackRepTime.h>
#include <HMatrixU.h>
#include <HMatrixV.h>
#include <HMatrixUV.h>
#include <HMatrixU6.h>
#include <HMatrixV6.h>
#include <HMatrixUV6.h>

#include <TEveElement.h>
#include <TEveGeoShape.h>
#include <TEveBox.h>
#include <TGeoManager.h>
#include <TGeoEltu.h>
#include <TMatrixDSymEigen.h>
#include <cassert>


namespace genfit {

PlanarMeasurement::PlanarMeasurement(int nDim)
  : AbsMeasurement(nDim), physicalPlane_(), planeId_(-1), stripV_(false)
{
  assert(nDim >= 1);
}

PlanarMeasurement::PlanarMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint), physicalPlane_(), planeId_(-1), stripV_(false)
{
  assert(rawHitCoords_.GetNrows() >= 1);
}


SharedPlanePtr PlanarMeasurement::constructPlane(const StateOnPlane&) const {
  if (!physicalPlane_) {
    Exception exc("PlanarMeasurement::constructPlane(): No plane has been set!", __LINE__,__FILE__);
    throw exc;
  }
  return physicalPlane_;
}


std::vector<MeasurementOnPlane*> PlanarMeasurement::constructMeasurementsOnPlane(const StateOnPlane& state) const {

  MeasurementOnPlane* mop = new MeasurementOnPlane(rawHitCoords_,
       rawHitCov_,
       state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mop);
  return retVal;
}


const AbsHMatrix* PlanarMeasurement::constructHMatrix(const AbsTrackRep* rep) const {

  if (dynamic_cast<const RKTrackRep*>(rep)
      || dynamic_cast<const RKTrackRepEnergy*>(rep) ) {

    switch(rawHitCoords_.GetNrows()) {
    case 1:
      if (stripV_)
        return new HMatrixV();
      return new HMatrixU();

    case 2:
      return new HMatrixUV();

    default:
      Exception exc("PlanarMeasurement default implementation can only handle 1D (strip) or 2D (pixel) measurements!", __LINE__,__FILE__);
      throw exc;
    }
  }

  if (dynamic_cast<const RKTrackRepTime*>(rep)) {
    switch(rawHitCoords_.GetNrows()) {
    case 1:
      if (stripV_)
        return new HMatrixV6();
      return new HMatrixU6();

    case 2:
      return new HMatrixUV6();

    default:
      Exception exc("PlanarMeasurement default implementation can only handle 1D (strip) or 2D (pixel) measurements!", __LINE__,__FILE__);
      throw exc;
    }
  }
    
  Exception exc("SpacepointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
  throw exc;
}

void PlanarMeasurement::drawDetector(TEveElementList *list) const
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

void PlanarMeasurement::drawMeasurement(TEveElementList *list, const MeasuredStateOnPlane& fittedState) const
{
  const TVector3& track_pos = fittedState.getPos();
  const TVector3& o = fittedState.getPlane()->getO();
  const TVector3& u = fittedState.getPlane()->getU();
  const TVector3& v = fittedState.getPlane()->getV();

  double hit_u = getRawHitCoords()(0);
  const TMatrixDSym& hit_cov = getRawHitCov();

  double errorScale_ = 1;
  double plane_size = 4;

  if(getDim() == 1) {
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
    return;
  }

  assert(getDim() == 2);

  double hit_v = getRawHitCoords()(1);

  // calculate eigenvalues to draw error-ellipse ----------------------------
  TMatrixDSymEigen eigen_values(getRawHitCov());
  TEveGeoShape* cov_shape = new TEveGeoShape("pixelCov");
  const TVectorD& ev = eigen_values.GetEigenValues();
  const TMatrixD& eVec = eigen_values.GetEigenVectors();
  double pseudo_res_0 = errorScale_*std::sqrt(ev(0));
  double pseudo_res_1 = errorScale_*std::sqrt(ev(1));
  // finished calcluating, got the values -----------------------------------
  
  // do autoscaling if necessary --------------------------------------------
  if(1 /*drawAutoScale_*/) {
    double min_cov = std::min(pseudo_res_0,pseudo_res_1);
    if(min_cov < 1e-5) {
      //std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
    } else {
      if(min_cov < 0.049) {
        double cor = 0.05 / min_cov;
        //std::cout << "Track " << i << ", Hit " << j << ": Pixel covariance too small, rescaling by " << cor;
        errorScale_ *= cor;
        pseudo_res_0 *= cor;
        pseudo_res_1 *= cor;
        //std::cout << " to " << errorScale_ << std::endl;
      }
    }
  }
  // finished autoscaling ---------------------------------------------------

  // calculate the semiaxis of the error ellipse ----------------------------
  cov_shape->SetShape(new TGeoEltu(pseudo_res_0, pseudo_res_1, 0.0105));
  TVector3 pix_pos = o + hit_u*u + hit_v*v;
  TVector3 u_semiaxis = (pix_pos + eVec(0,0)*u + eVec(1,0)*v)-pix_pos;
  TVector3 v_semiaxis = (pix_pos + eVec(0,1)*u + eVec(1,1)*v)-pix_pos;
  TVector3 norm = u.Cross(v);
  // finished calculating ---------------------------------------------------

  // rotate and translate everything correctly ------------------------------
  TGeoRotation det_rot("det_rot", (u_semiaxis.Theta()*180)/TMath::Pi(), (u_semiaxis.Phi()*180)/TMath::Pi(),
                       (v_semiaxis.Theta()*180)/TMath::Pi(), (v_semiaxis.Phi()*180)/TMath::Pi(),
                       (norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi());
  TGeoCombiTrans det_trans(pix_pos(0),pix_pos(1),pix_pos(2), &det_rot);
  cov_shape->SetTransMatrix(det_trans);
  // finished rotating and translating --------------------------------------

  cov_shape->SetMainColor(kYellow);
  cov_shape->SetMainTransparency(0);

  list->AddElement(cov_shape);
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
