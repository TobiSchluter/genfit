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

#ifndef  _MeasurementOnPlane_h
#define  _MeasurementOnPlane_h

#include <TMatrixD.h>

#include "MeasuredStateOnPlane.h"



namespace genfit {

class MeasurementOnPlane : public MeasuredStateOnPlane {

 public:

  MeasurementOnPlane();
  MeasurementOnPlane(const TVectorD& state, const TMatrixDSym& cov, DetPlane* plane, AbsTrackRep* rep, const TMatrixD& hMatrix, double weight = 1.);

  const TMatrixD& getHMatrix() const {return hMatrix_;}
  double getWeight() const {return weight_;}

  void setHMatrix(const TMatrixD& hMatrix) {hMatrix_ = hMatrix;}
  void setWeight(double weight) {weight_ = weight;}

 protected:

  TMatrixD hMatrix_; // projection matrix
  double weight_;


  ClassDef(MeasurementOnPlane,1)

};

} /* End of namespace   */

#endif //  _MeasurementOnPlane_h
