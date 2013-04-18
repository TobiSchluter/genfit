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

namespace genfit {

ReferenceStateOnPlane::ReferenceStateOnPlane() :
  StateOnPlane(),
  forwardSegmentLength_(0),
  backwardSegmentLength_(0),
  forwardTransportMatrix_(),
  backwardTransportMatrix_(),
  forwardNoiseMatrix_(),
  backwardNoiseMatrix_()
{
  ;
}

ReferenceStateOnPlane::ReferenceStateOnPlane(const TVectorD& state,
    sharedPlanePtr plane,
    AbsTrackRep* rep,
    double forwardSegmentLength,
    double backwardSegmentLength,
    const TMatrixD& forwardTransportMatrix,
    const TMatrixD& backwardTransportMatrix,
    const TMatrixDSym& forwardNoiseMatrix,
    const TMatrixDSym& backwardNoiseMatrix) :
  StateOnPlane(state, plane, rep),
  forwardSegmentLength_(forwardSegmentLength),
  backwardSegmentLength_(backwardSegmentLength),
  forwardTransportMatrix_(forwardTransportMatrix),
  backwardTransportMatrix_(backwardTransportMatrix),
  forwardNoiseMatrix_(forwardNoiseMatrix),
  backwardNoiseMatrix_(backwardNoiseMatrix)
{
  ;
}

ReferenceStateOnPlane::ReferenceStateOnPlane(const StateOnPlane& state,
    double forwardSegmentLength,
    double backwardSegmentLength,
    const TMatrixD& forwardTransportMatrix,
    const TMatrixD& backwardTransportMatrix,
    const TMatrixDSym& forwardNoiseMatrix,
    const TMatrixDSym& backwardNoiseMatrix) :
  StateOnPlane(state),
  forwardSegmentLength_(forwardSegmentLength),
  backwardSegmentLength_(backwardSegmentLength),
  forwardTransportMatrix_(forwardTransportMatrix),
  backwardTransportMatrix_(backwardTransportMatrix),
  forwardNoiseMatrix_(forwardNoiseMatrix),
  backwardNoiseMatrix_(backwardNoiseMatrix)
{
  ;
}

} /* End of namespace genfit */
