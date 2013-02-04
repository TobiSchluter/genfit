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

#include "GFStateOnPlane5D.h"

#include <assert.h>
#include <iostream>
#include <limits>


GFStateOnPlane5D::GFStateOnPlane5D()
  : GFAbsStateOnPlane(DIM),
    fCharge(0),
    fSpu(true)
{
  ;
}

GFStateOnPlane5D::GFStateOnPlane5D(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, double charge, bool spu)
  : GFAbsStateOnPlane(plane, state),
    fCharge(charge),
    fSpu(spu)
{
  assert(state.GetNrows() == DIM);
}


void GFStateOnPlane5D::setData(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, double charge, bool spu){
  GFAbsStateOnPlane::setData(plane, state);
  fCharge = charge;
  fSpu = spu;
}

TVector3 GFStateOnPlane5D::getPos() const {

  const TVector3& U(fPlane->getU());
  const TVector3& V(fPlane->getV());
  const TVector3& O(fPlane->getO());

  return TVector3(O.X() + fState(3)*U.X() + fState(4)*V.X(),  // x
                  O.Y() + fState(3)*U.Y() + fState(4)*V.Y(),  // y
                  O.Z() + fState(3)*U.Z() + fState(4)*V.Z()); // z

}

TVector3 GFStateOnPlane5D::getMom() const{

  const TVector3& U(fPlane->getU());
  const TVector3& V(fPlane->getV());
  TVector3 W(fPlane->getNormal());

  TVector3 mom(W.X() + fState(1)*U.X() + fState(2)*V.X(),  // a_x
               W.Y() + fState(1)*U.Y() + fState(2)*V.Y(),  // a_y
               W.Z() + fState(1)*U.Z() + fState(2)*V.Z()); // a_z

  mom.SetMag(fCharge/fState(0));
  if (!fSpu) mom *= -1.;

  return mom;

}

void GFStateOnPlane5D::getPosMom(TVector3& pos, TVector3& mom) const {

  const TVector3& U(fPlane->getU());
  const TVector3& V(fPlane->getV());
  const TVector3& O(fPlane->getO());
  TVector3 W(fPlane->getNormal());

  pos.SetXYZ(O.X() + fState(3)*U.X() + fState(4)*V.X(),  // x
             O.Y() + fState(3)*U.Y() + fState(4)*V.Y(),  // y
             O.Z() + fState(3)*U.Z() + fState(4)*V.Z()); // z

  mom.SetXYZ(W.X() + fState(1)*U.X() + fState(2)*V.X(),  // a_x
             W.Y() + fState(1)*U.Y() + fState(2)*V.Y(),  // a_y
             W.Z() + fState(1)*U.Z() + fState(2)*V.Z()); // a_z

  mom.SetMag(fCharge/fState(0));
  if (!fSpu) mom *= -1.;

}

