/* Copyright 2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Johannes Rauch, Tobias Schlüter

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

#include "HMatrixV6.h"
#include <cassert>
#include <alloca.h>
#include <iostream>

namespace genfit {


// 0, 0, 0, 0, 1, 0

const TMatrixD& HMatrixV6::getMatrix() const {
  static const double HMatrixContent[6] = {0, 0, 0, 0, 1, 0};

  static const TMatrixD HMatrix(1,6, HMatrixContent);

  return HMatrix;
}


TVectorD HMatrixV6::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 6);

  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = v(4); // v

  return TVectorD(1, retValArray);
}


TMatrixD HMatrixV6::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 6);

  double* retValArray =(double *)alloca(sizeof(double) * 6);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<6; ++i) {
    retValArray[i] = MatArray[i*6 + 4];
  }

  return TMatrixD(6,1, retValArray);
}


TMatrixD HMatrixV6::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 6);

  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows());
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i] = MatArray[i*6 + 4];
  }

  return TMatrixD(M.GetNrows(),1, retValArray);
}


void HMatrixV6::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 6);

  M(0,0) = M(4,4);

  M.ResizeTo(1,1);
}


void HMatrixV6::Print(const Option_t*) const {
  std::cout << "V" << std::endl;
}

} /* End of namespace genfit */
