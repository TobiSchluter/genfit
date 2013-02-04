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

#include "GFAbsStateOnPlane.h"

#include <assert.h>
#include <iostream>
#include <limits>


GFAbsStateOnPlane::GFAbsStateOnPlane(unsigned int dim)
  : fPlane(new GFDetPlane()),
    fState(dim)
{
  ;
}

GFAbsStateOnPlane::GFAbsStateOnPlane(SHARED_PTR(GFDetPlane) plane, const TVectorD& state)
  : fPlane(plane),
    fState(state)
{
  ;
}


void GFAbsStateOnPlane::setData(SHARED_PTR(GFDetPlane) plane, const TVectorD& state) {
  fPlane = plane;
  setState(state);
}

void GFAbsStateOnPlane::setState(const TVectorD& state) {
  assert(state.GetNrows() == fState.GetNrows());
  fState = state;
}


void GFAbsStateOnPlane::Print(const Option_t*) const {
  std::cout << "GFAbsStateOnPlane defined in plane:"; fPlane->Print();
  std::cout << "state vector: "; fState.Print();
}
