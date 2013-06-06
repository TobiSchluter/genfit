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

/** @addtogroup genfit
 * @{
 */

#ifndef genfit_ReferenceStateOnPlane_h
#define genfit_ReferenceStateOnPlane_h

#include "StateOnPlane.h"


namespace genfit {


  /** 
   *  Transport matrices describe transport TO that plane.
   */
class ReferenceStateOnPlane : public StateOnPlane {

 public:

  ReferenceStateOnPlane();
  ReferenceStateOnPlane(const TVectorD& state,
      SharedPlanePtr plane,
      AbsTrackRep* rep,
      double forwardSegmentLength,
      double backwardSegmentLength,
      const TMatrixD& forwardTransportMatrix,
      const TMatrixD& backwardTransportMatrix,
      const TMatrixDSym& forwardNoiseMatrix,
      const TMatrixDSym& backwardNoiseMatrix);
  ReferenceStateOnPlane(const StateOnPlane& state,
      double forwardSegmentLength,
      double backwardSegmentLength,
      const TMatrixD& forwardTransportMatrix,
      const TMatrixD& backwardTransportMatrix,
      const TMatrixDSym& forwardNoiseMatrix,
      const TMatrixDSym& backwardNoiseMatrix);

  double getForwardSegmentLength() const {return forwardSegmentLength_;}
  double getBackwardSegmentLength() const {return backwardSegmentLength_;}
  const TMatrixD& getForwardTransportMatrix() const {return forwardTransportMatrix_;}
  const TMatrixD& getBackwardTransportMatrix() const {return backwardTransportMatrix_;}
  const TMatrixDSym& getForwardNoiseMatrix() const {return forwardNoiseMatrix_;}
  const TMatrixDSym& getBackwardNoiseMatrix() const {return backwardNoiseMatrix_;}

  virtual void Print(Option_t* option = "") const override;

 protected:

  double forwardSegmentLength_; // Segment length from previous referenceState
  double backwardSegmentLength_; //  Segment length from next referenceState
  TMatrixD forwardTransportMatrix_; // transport matrix from previous referenceState
  TMatrixD backwardTransportMatrix_; // transport matrix from next referenceState
  TMatrixDSym forwardNoiseMatrix_; // noise matrix for transport from previous referenceState
  TMatrixDSym backwardNoiseMatrix_; // noise matrix for transport from next referenceState


  //ClassDef(ReferenceStateOnPlane,0)

};

} /* End of namespace genfit */

#endif // genfit_ReferenceStateOnPlane_h
