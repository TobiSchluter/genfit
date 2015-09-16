/* Copyright 2015, Ludwig-Maximilians-Universität München,
   Authors: Tobias Schlüter

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

#include "HMatrixUT6.h"
#include <cassert>
#include <alloca.h>
#include <iostream>


namespace genfit {


// 0, 0, 0, 1, 0, 0
// 0, 0, 0, 0, 0, 1

const TMatrixD& HMatrixUT6::getMatrix() const {
  static const double HMatrixContent[2*6] = {0, 0, 0, 1, 0, 0,
                                             0, 0, 0, 0, 0, 1};

  static const TMatrixD HMatrix(2,6, HMatrixContent);

  return HMatrix;
}


TVectorD HMatrixUT6::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 6);

  double* retValArray =(double *)alloca(sizeof(double) * 2);

  retValArray[0] = v(3); // u
  retValArray[1] = v(5); // t

  return TVectorD(2, retValArray);
}


TMatrixD HMatrixUT6::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 6);

  double* retValArray =(double *)alloca(sizeof(double) * 2 * 6);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<6; ++i) {
    retValArray[i*2    ] = MatArray[i*6 + 3];
    retValArray[i*2 + 1] = MatArray[i*6 + 5];
  }

  return TMatrixD(6,2, retValArray);
}


TMatrixD HMatrixUT6::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 6);

  double* retValArray =(double *)alloca(sizeof(double) * 2 * M.GetNrows());
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i*2    ] = MatArray[i*6 + 3];
    retValArray[i*2 + 1] = MatArray[i*6 + 5];
  }

  return TMatrixD(M.GetNrows(),2, retValArray);
}


void HMatrixUT6::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 6);

  M(0,0) = M(3,3);
  M(1,0) = M(0,1) = M(3,5);
  M(1,1) = M(5,5);

  M.ResizeTo(2,2);
}


void HMatrixUT6::Print(const Option_t*) const {
  std::cout << "UT6" << std::endl;
}


} /* End of namespace genfit */
