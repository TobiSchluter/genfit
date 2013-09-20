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

  double* retValArray =(double *)alloca(sizeof(double) * 2);
  const double* VecArray = v.GetMatrixArray();

  retValArray[0] = VecArray[3]; // u
  retValArray[1] = VecArray[4]; // v

  return TVectorD(2, retValArray);
}


TMatrixD HMatrixUV::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 5*2);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i*2] = MatArray[i*5 + 3];
    retValArray[i*2 + 1] = MatArray[i*5 + 4];
  }

  return TMatrixD(5,2, retValArray);
}


TMatrixD HMatrixUV::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows()*2);
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i*2] = MatArray[i*5 + 3];
    retValArray[i*2 + 1] = MatArray[i*5 + 4];
  }

  return TMatrixD(M.GetNrows(),2, retValArray);
}


void HMatrixUV::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  double* MatArray = M.GetMatrixArray();

  for (unsigned int i = 0; i < 2; ++i) {
    MatArray[i] = MatArray[18+i];
    MatArray[5+i] = MatArray[23+i];
  }

  M.ResizeTo(2,2);
}


} /* End of namespace genfit */
