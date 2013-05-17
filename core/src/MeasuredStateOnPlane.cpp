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

#include "MeasuredStateOnPlane.h"
#include "AbsTrackRep.h"

#include <assert.h>
#include <iostream>

namespace genfit {

MeasuredStateOnPlane::MeasuredStateOnPlane(const AbsTrackRep* rep) :
  StateOnPlane(rep), cov_(0,0)
{
  if (rep != nullptr) {
    cov_.ResizeTo(rep->getDim(), rep->getDim());
  }
}

MeasuredStateOnPlane::MeasuredStateOnPlane(const TVectorD& state, const TMatrixDSym& cov, SharedPlanePtr plane, const AbsTrackRep* rep) :
  StateOnPlane(state, plane, rep), cov_(cov)
{
  assert(rep != nullptr);
  //assert(cov_.GetNcols() == (signed)rep->getDim());
}

MeasuredStateOnPlane::MeasuredStateOnPlane(const StateOnPlane& state, const TMatrixDSym& cov) :
  StateOnPlane(state), cov_(cov)
{
  //assert(cov_.GetNcols() == (signed)getRep()->getDim());
}

void MeasuredStateOnPlane::Print(Option_t* option) const {
  std::cout << "genfit::MeasuredStateOnPlane ";
  std::cout << " state vector: "; state_.Print();
  std::cout << " covariance matrix: "; cov_.Print();
  if (sharedPlane_ != nullptr) {
    std::cout << " defined in plane "; sharedPlane_->Print();
    TVector3 pos, mom;
    TMatrixDSym cov(6,6);
    getRep()->getPosMomCov(this, pos, mom, cov);
    std::cout << " 3D position: "; pos.Print();
    std::cout << " 3D momentum: "; mom.Print();
    std::cout << " 6D covariance: "; cov.Print();
  }
}

} /* End of namespace genfit */
