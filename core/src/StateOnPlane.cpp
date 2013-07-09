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

#include <cassert>
#include <iostream>

namespace genfit {

StateOnPlane::StateOnPlane(const AbsTrackRep* rep) :
  state_(0), auxInfo_(0), sharedPlane_(), rep_(rep)
{
  if (rep != NULL) {
    state_.ResizeTo(rep->getDim());
  }
}

StateOnPlane::StateOnPlane(const TVectorD& state, SharedPlanePtr plane, const AbsTrackRep* rep) :
  state_(state), sharedPlane_(plane), rep_(rep)
{
  assert(rep != NULL);
  //assert(state_.GetNrows() == (signed)rep->getDim());
}

StateOnPlane& StateOnPlane::operator= (const StateOnPlane& other) {
  state_.ResizeTo(other.state_);
  state_ = other.state_;

  auxInfo_.ResizeTo(other.auxInfo_);
  auxInfo_ = other.auxInfo_;

  sharedPlane_ = other.sharedPlane_;

  rep_ = other.rep_;

  return *this;
}


void StateOnPlane::Print(Option_t* option) const {
  std::cout << "genfit::StateOnPlane ";
  std::cout << " state vector: "; state_.Print();
  if (sharedPlane_ != NULL) {
    std::cout << " defined in plane "; sharedPlane_->Print();
    TVector3 pos, mom;
    getRep()->getPosMom(this, pos, mom);
    std::cout << " 3D position: "; pos.Print();
    std::cout << " 3D momentum: "; mom.Print();
  }
}

} /* End of namespace genfit */
