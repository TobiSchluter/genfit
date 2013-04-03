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
      const DetPlane* plane,
      const AbsTrackRep* rep,
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


 protected:

  double forwardSegmentLength_;
  double backwardSegmentLength_;
  TMatrixD forwardTransportMatrix_;
  TMatrixD backwardTransportMatrix_;
  TMatrixDSym forwardNoiseMatrix_;
  TMatrixDSym backwardNoiseMatrix_;


  ClassDef(ReferenceStateOnPlane,0)

};

} /* End of namespace genfit */

#endif // genfit_ReferenceStateOnPlane_h
