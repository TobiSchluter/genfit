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

namespace genfit {

StateOnPlane::StateOnPlane() :
  state_(0), sharedPlane_(nullptr), rep_(nullptr)
{
  ;
}

StateOnPlane::StateOnPlane(const TVectorD& state, const DetPlane* plane, const AbsTrackRep* rep) :
  state_(state), sharedPlane_(plane), rep_(rep)
{
  ;
}


StateOnPlane::~StateOnPlane() {
  if (sharedPlane_ != nullptr)
    delete sharedPlane_;
}


void
StateOnPlane::setStatePlane(const TVectorD& state, const DetPlane* plane) {
  state_ = state;
  if (sharedPlane_ != nullptr)
    delete sharedPlane_;
  sharedPlane_ = plane;
}

} /* End of namespace genfit */
