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

#include "GFAbsMeasuredStateOnPlane.h"

#include <assert.h>
#include <iostream>
#include <limits>


GFAbsMeasuredStateOnPlane::GFAbsMeasuredStateOnPlane(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, const TMatrixDSym& cov)
  : GFAbsStateOnPlane(plane, state), fCov(cov)
{
  checkDim();
}

GFAbsMeasuredStateOnPlane::GFAbsMeasuredStateOnPlane(const GFAbsStateOnPlane& state, const TMatrixDSym& cov)
  : GFAbsStateOnPlane(state), fCov(cov)
{
  checkDim();
}


void GFAbsMeasuredStateOnPlane::setData(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, const TMatrixDSym& cov){
  GFAbsStateOnPlane::setData(plane, state);
  setCov(cov);
}

void GFAbsMeasuredStateOnPlane::setStateCov(const TVectorD& state, const TMatrixDSym& cov){
  setState(state);
  setCov(cov);
}

void GFAbsMeasuredStateOnPlane::setCov(const TMatrixDSym& cov){
  fCov = cov;
  checkDim();
}


void GFAbsMeasuredStateOnPlane::Print(const Option_t*) const {
  GFAbsStateOnPlane::Print();
  std::cout << "covariance matrix: "; fCov.Print();
}


void GFAbsMeasuredStateOnPlane::checkDim(){
  assert(fState.GetNrows() == fCov.GetNrows());
}
