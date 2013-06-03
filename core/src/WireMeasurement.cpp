/* Copyright 2013, Ludwig-Maximilians Universitaet Muenchen,
   Authors: Tobias Schlüter

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

#include <math.h>
#include <TMatrixD.h>

#include "AbsTrackRep.h"
#include "RKTrackRep.h"
#include "MeasurementOnPlane.h"
#include "Exception.h"

#include "WireMeasurement.h"

namespace genfit {


WireMeasurement::WireMeasurement(int nDim)
  : AbsMeasurement(nDim), maxDistance_(1.E50), leftRight_(0)
{
  ;
}

WireMeasurement::WireMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint), maxDistance_(1.E50), leftRight_(0)
{
  ;
}

SharedPlanePtr WireMeasurement::constructPlane(const StateOnPlane* state) const {

  // copy state. Neglect covariance.
  StateOnPlane st(*state);

  TVector3 wire1(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));
  TVector3 wire2(rawHitCoords_(3), rawHitCoords_(4), rawHitCoords_(5));

  //std::cout << " wire1(" << rawHitCoords_(0) << ", " << rawHitCoords_(1) << ", " << rawHitCoords_(2) << ")" << std::endl;
  //std::cout << " wire2(" << rawHitCoords_(3) << ", " << rawHitCoords_(4) << ", " << rawHitCoords_(5) << ")" << std::endl;

  // unit vector along the wire (V)
  TVector3 wireDirection = wire2 - wire1; 
  wireDirection.SetMag(1.);

  //std::cout << " wireDirection(" << wireDirection.X() << ", " << wireDirection.Y() << ", " << wireDirection.Z() << ")" << std::endl;

  // point of closest approach
  const AbsTrackRep* rep = state->getRep();
  rep->extrapolateToLine(&st, wire1, wireDirection);
  const TVector3& poca = rep->getPos(&st);
  TVector3 dirInPoca = rep->getMom(&st);
  dirInPoca.SetMag(1.);
  const TVector3& pocaOnWire = wire1 + wireDirection.Dot(poca - wire1)*wireDirection;


  // check distance of poca to wire
  if((poca - pocaOnWire).Mag() > maxDistance_) {
    Exception exc("GFAbsWireHit::detPlane(): distance poca-wire > maxdistance", __LINE__,__FILE__);
    throw exc;    
  }

 
  // check if direction is parallel to wire
  if (fabs(wireDirection.Angle(dirInPoca)) < 0.01){
    Exception exc("GFAbsWireHit::detPlane(): Cannot construct detector plane, direction is parallel to wire", __LINE__,__FILE__);
    throw exc;
  }
  
  // construct orthogonal vector
  TVector3 U = dirInPoca.Cross(wireDirection);
  // U.SetMag(1.); automatically assured

  if ((poca - pocaOnWire)*U < 0)
    U *= -1.;

  // check left/right ambiguity
  if (leftRight_ == 0){ // auto select
    if ((poca - pocaOnWire)*U < 0) U *= -1.;
  }
  else if (leftRight_ < 0) U *= -1.;

  return SharedPlanePtr(new DetPlane(wire1, U, wireDirection, 0));
}


MeasurementOnPlane WireMeasurement::constructMeasurementOnPlane(const AbsTrackRep* rep, const SharedPlanePtr plane) const
{
  double m = rawHitCoords_(6);
  double V = rawHitCov_(6,6);
  MeasurementOnPlane mop(TVectorD(1, &m),
			 TMatrixDSym(1, &V),
			 plane, rep, getHMatrix(rep));
  return mop;
}


const TMatrixD& WireMeasurement::getHMatrix(const AbsTrackRep* rep) const {
  if (dynamic_cast<const RKTrackRep*>(rep) != nullptr) {
    static const double HMatrixContent[5] = {0, 0, 0, 1, 0};
    static const TMatrixT<double> HMatrix(1,5, HMatrixContent);

    return HMatrix;
  }
  else {
    Exception exc("WireMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }
}


void WireMeasurement::setLeftRightResolution(int lr){
  if (lr==0) leftRight_ = 0;
  else if (lr<0) leftRight_ = -1;
  else leftRight_ = 1;
}


} /* End of namespace genfit */
