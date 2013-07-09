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

#include <iostream>

#include "MeasurementOnPlane.h"

namespace genfit {

MeasurementOnPlane::MeasurementOnPlane(const AbsTrackRep* rep) :
  MeasuredStateOnPlane(rep), hMatrix_(0,0), weight_(0)
{
  ;
}

MeasurementOnPlane::MeasurementOnPlane(const TVectorD& state, const TMatrixDSym& cov, SharedPlanePtr plane, const AbsTrackRep* rep, const TMatrixD& hMatrix, double weight) :
  MeasuredStateOnPlane(state, cov, plane, rep), hMatrix_(hMatrix), weight_(weight)
{
  ;
}

void MeasurementOnPlane::Print(Option_t* option) const
{
  std::cout << "genfit::MeasurementOnPlane, weight = " << weight_ << "\n";
  std::cout << " state vector: "; state_.Print();
  std::cout << " covariance matrix: "; cov_.Print();
  if (sharedPlane_ != NULL)
    std::cout << " defined in plane "; sharedPlane_->Print();
  std::cout << " hMatrix: "; hMatrix_.Print();

}

} /* End of namespace   */
