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

#include "StateOnPlane.h"
#include "AbsTrackRep.h"

#include <assert.h>
#include <iostream>

namespace genfit {

StateOnPlane::StateOnPlane(AbsTrackRep* rep) :
  state_(0), auxInfo_(0), sharedPlane_(), rep_(rep)
{
  if (rep != nullptr) {
    state_.ResizeTo(rep->getDim());
  }
}

StateOnPlane::StateOnPlane(const TVectorD& state, SharedPlanePtr plane, AbsTrackRep* rep) :
  state_(state), sharedPlane_(plane), rep_(rep)
{
  assert(rep != nullptr);
  assert(state_.GetNrows() == rep->getDim());
}

void StateOnPlane::Print(Option_t* option) const {
  std::cout << "genfit::StateOnPlane ";
  std::cout << " state vector: "; state_.Print();
  if (sharedPlane_ != nullptr) {
    std::cout << " defined in plane "; sharedPlane_->Print();
    TVector3 pos, mom;
    getRep()->getPosMom(this, pos, mom);
    std::cout << " 3D position: "; pos.Print();
    std::cout << " 3D momentum: "; mom.Print();
  }
}

} /* End of namespace genfit */
