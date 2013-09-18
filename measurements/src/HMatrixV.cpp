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

#include <HMatrixV.h>


namespace genfit {


// 0, 0, 0, 0, 1

TVectorD HMatrixV::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 5);

  TVectorD retVal(1);
  retVal(0) = v(4); // u

  return retVal;
}


TMatrixD HMatrixV::MHt(const TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  TMatrixD retVal(5,1);

  retVal(0,0) = M(0,4);
  retVal(1,0) = M(1,4);
  retVal(2,0) = M(2,4);
  retVal(3,0) = M(3,4);
  retVal(4,0) = M(4,4);

  return retVal;
}


TMatrixD HMatrixV::MHt(const TMatrixD& M) const {
  assert (M.GetNrows() == 5);

  TMatrixD retVal(5,1);

  retVal(0,0) = M(0,4);
  retVal(1,0) = M(1,4);
  retVal(2,0) = M(2,4);
  retVal(3,0) = M(3,4);
  retVal(4,0) = M(4,4);

  return retVal;
}


void HMatrixV::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  M(0,0) = M(4,4);

  M.ResizeTo(1,1);
}


} /* End of namespace genfit */
