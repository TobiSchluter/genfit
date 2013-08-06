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

#include "ReferenceStateOnPlane.h"

#include <iostream>


namespace genfit {

ReferenceStateOnPlane::ReferenceStateOnPlane() :
  StateOnPlane(),
  forwardSegmentLength_(0),
  backwardSegmentLength_(0),
  forwardTransportMatrix_(),
  backwardTransportMatrix_(),
  forwardNoiseMatrix_(),
  backwardNoiseMatrix_(),
  forwardDeltaState_(),
  backwardDeltaState_()
{
  ;
}

ReferenceStateOnPlane::ReferenceStateOnPlane(const TVectorD& state,
    SharedPlanePtr plane,
    const AbsTrackRep* rep) :
  StateOnPlane(state, plane, rep),
  forwardSegmentLength_(0),
  backwardSegmentLength_(0),
  forwardTransportMatrix_(rep->getDim(), rep->getDim()),
  backwardTransportMatrix_(rep->getDim(), rep->getDim()),
  forwardNoiseMatrix_(rep->getDim()),
  backwardNoiseMatrix_(rep->getDim()),
  forwardDeltaState_(rep->getDim()),
  backwardDeltaState_(rep->getDim())
{
  ;
}

ReferenceStateOnPlane::ReferenceStateOnPlane(const TVectorD& state,
    SharedPlanePtr plane,
    const AbsTrackRep* rep,
    const TVectorD& auxInfo) :
  StateOnPlane(state, plane, rep, auxInfo),
  forwardSegmentLength_(0),
  backwardSegmentLength_(0),
  forwardTransportMatrix_(rep->getDim(), rep->getDim()),
  backwardTransportMatrix_(rep->getDim(), rep->getDim()),
  forwardNoiseMatrix_(rep->getDim()),
  backwardNoiseMatrix_(rep->getDim()),
  forwardDeltaState_(rep->getDim()),
  backwardDeltaState_(rep->getDim())
{
  ;
}


ReferenceStateOnPlane::ReferenceStateOnPlane(const StateOnPlane& state) :
  StateOnPlane(state),
  forwardSegmentLength_(0),
  backwardSegmentLength_(0),
  forwardTransportMatrix_(state.getRep()->getDim(), state.getRep()->getDim()),
  backwardTransportMatrix_(state.getRep()->getDim(), state.getRep()->getDim()),
  forwardNoiseMatrix_(state.getRep()->getDim()),
  backwardNoiseMatrix_(state.getRep()->getDim()),
  forwardDeltaState_(state.getRep()->getDim()),
  backwardDeltaState_(state.getRep()->getDim())
{
  std::cout << "should never come here" << std::endl;
  exit(0);
}


void ReferenceStateOnPlane::Print(Option_t* option) const {
  StateOnPlane::Print();

  std::cout << " forwardSegmentLength: " << forwardSegmentLength_ << "\n";
  std::cout << " forwardTransportMatrix: "; forwardTransportMatrix_.Print();
  std::cout << " forwardNoiseMatrix: "; forwardNoiseMatrix_.Print();
  std::cout << " forwardDeltaState: "; forwardDeltaState_.Print();

  std::cout << " backwardSegmentLength_: " << backwardSegmentLength_ << "\n";
  std::cout << " backwardTransportMatrix: "; backwardTransportMatrix_.Print();
  std::cout << " backwardNoiseMatrix: "; backwardNoiseMatrix_.Print();
  std::cout << " backwardDeltaState: "; backwardDeltaState_.Print();

}


} /* End of namespace genfit */
