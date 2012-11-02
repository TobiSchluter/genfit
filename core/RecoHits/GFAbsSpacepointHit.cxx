/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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
#include "GFAbsSpacepointHit.h"
#include "assert.h"
#include "TMath.h"


void
GFAbsSpacepointHit::getMeasurement(const GFAbsTrackRep* rep,
                                const GFDetPlane& pl,
                                const TMatrixT<double>& statePred,
                                const TMatrixT<double>& covPred,
                                TMatrixT<double>& m,
                                TMatrixT<double>& V) {

  static_cast<void>(rep);
  static_cast<void>(statePred);
  static_cast<void>(covPred);

  TVector3 o(pl.getO());
  TVector3 u(pl.getU());
  TVector3 v(pl.getV());

  // m
  m.ResizeTo(2,1);

  TMatrixT<double> D(3,1);
  D(0,0) = o.X();
  D(1,0) = o.Y();
  D(2,0) = o.Z();

  D *= -1.;
  D += fHitCoord;
  //now the vector D points from the origin of the plane to the hit point

  m(0,0) = D(0,0) * u.X() + D(1,0) * u.Y() + D(2,0) * u.Z();
  m(1,0) = D(0,0) * v.X() + D(1,0) * v.Y() + D(2,0) * v.Z();


  // V
  V.ResizeTo(2,2);

  TMatrixT<double> jac(3,2);
  
  // jac = dF_i/dx_j = s_unitvec * t_untivec, with s=u,v and t=x,y,z
  jac(0,0) = u.X();
  jac(1,0) = u.Y();
  jac(2,0) = u.Z();
  jac(0,1) = v.X();
  jac(1,1) = v.Y();
  jac(2,1) = v.Z();

  TMatrixT<double> jac_orig = jac;
  TMatrixT<double> jac_t = jac.T();

  V = jac_t * (fHitCov * jac_orig); // TODO use similarity

}


const GFDetPlane&
GFAbsSpacepointHit::getDetPlane(GFAbsTrackRep* rep)
{
  TVector3 point(fHitCoord(0,0), fHitCoord(1,0), fHitCoord(2,0));

  TVector3 poca, dirInPoca;
  rep->extrapolateToPoint(point, poca, dirInPoca);

  fPlane.setO(point);
  fPlane.setNormal(dirInPoca);

  return fPlane;
}


ClassImp(GFAbsSpacepointHit)
