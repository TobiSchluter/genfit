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

#include <HMatrixUV.h>


namespace genfit {


// 0, 0, 0, 1, 0
// 0, 0, 0, 0, 1

TVectorD HMatrixUV::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 5);

  TVectorD retVal(2);
  retVal(0) = v(3); // u
  retVal(1) = v(4); // v

  return retVal;
}


TMatrixD HMatrixUV::MHt(const TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  TMatrixD retVal(5,2);

  retVal(0,0) = M(0,3);
  retVal(1,0) = M(1,3);
  retVal(2,0) = M(2,3);
  retVal(3,0) = M(3,3);
  retVal(4,0) = M(4,3);

  retVal(0,1) = M(0,4);
  retVal(1,1) = M(1,4);
  retVal(2,1) = M(2,4);
  retVal(3,1) = M(3,4);
  retVal(4,1) = M(4,4);

  return retVal;
}


TMatrixD HMatrixUV::MHt(const TMatrixD& M) const {
  assert (M.GetNrows() == 5);

  TMatrixD retVal(5,2);

  retVal(0,0) = M(0,3);
  retVal(1,0) = M(1,3);
  retVal(2,0) = M(2,3);
  retVal(3,0) = M(3,3);
  retVal(4,0) = M(4,3);

  retVal(0,1) = M(0,4);
  retVal(1,1) = M(1,4);
  retVal(2,1) = M(2,4);
  retVal(3,1) = M(3,4);
  retVal(4,1) = M(4,4);

  return retVal;
}


void HMatrixUV::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  M(0,0) = M(3,3);
  M(0,1) = M(3,4);
  M(1,0) = M(4,3);
  M(1,1) = M(4,4);

  M.ResizeTo(2,2);
}


} /* End of namespace genfit */
