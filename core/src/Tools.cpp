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

#include "Tools.h"

#include <cmath>
#include <memory>
#include <typeinfo>

#include <TDecompChol.h>
#include <TMath.h>

#include "Exception.h"

namespace genfit {

void tools::invertMatrix(const TMatrixDSym& mat, TMatrixDSym& inv, double* determinant){
  inv.ResizeTo(mat);

  // check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
  if (!(mat<1.E100) || !(mat>-1.E100)){
    Exception e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1e100)",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  // do the trivial inversions for 1x1 and 2x2 matrices manually
  if (mat.GetNrows() == 1){
    if (determinant != NULL) *determinant = mat(0,0);
    inv(0,0) = 1./mat(0,0);
    return;
  }

  if (mat.GetNrows() == 2){
    double det = mat(0,0)*mat(1,1) - mat(1,0)*mat(1,0);
    if (determinant != NULL) *determinant = det;
    if(fabs(det) < 1E-50){
      Exception e("cannot invert matrix GFTools::invertMatrix(), determinant = 0",
          __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    det = 1./det;
    inv(0,0) =             det * mat(1,1);
    inv(0,1) = inv(1,0) = -det * mat(1,0);
    inv(1,1) =             det * mat(0,0);
    return;
  }

  // else use TDecompChol
  bool status = 0;
  TDecompChol invertAlgo(mat, 1E-50);

  status = invertAlgo.Invert(inv);
  if(status == 0){
    Exception e("cannot invert matrix GFTools::invertMatrix(), status = 0",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  if (determinant != NULL) {
    double d1, d2;
    invertAlgo.Det(d1, d2);
    *determinant = ldexp(d1, d2);
  }
}

void tools::invertMatrix(TMatrixDSym& mat, double* determinant){
  // check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
  if (!(mat<1.E100) || !(mat>-1.E100)){
    Exception e("cannot invert matrix GFTools::invertMatrix(), entries too big (>1e100)",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  // do the trivial inversions for 1x1 and 2x2 matrices manually
  if (mat.GetNrows() == 1){
    if (determinant != NULL) *determinant = mat(0,0);
    mat(0,0) = 1./mat(0,0);
    return;
  }
  
  if (mat.GetNrows() == 2){
    double *arr = mat.GetMatrixArray();
    double det = arr[0]*arr[3] - arr[1]*arr[1];
    if (determinant != NULL) *determinant = det;
    if(fabs(det) < 1E-50){
      Exception e("cannot invert matrix GFTools::invertMatrix(), determinant = 0",
          __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    det = 1./det;
    double temp[3];
    temp[0] =  det * arr[3];
    temp[1] = -det * arr[1];
    temp[2] =  det * arr[0];
    //double *arr = mat.GetMatrixArray();
    arr[0] = temp[0];
    arr[1] = arr[2] = temp[1];
    arr[3] = temp[2];
    return;
  }

  // else use TDecompChol
  bool status = 0;
  TDecompChol invertAlgo(mat, 1E-50);

  status = invertAlgo.Invert(mat);
  if(status == 0){
    Exception e("cannot invert matrix GFTools::invertMatrix(), status = 0",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  if (determinant != NULL) {
    double d1, d2;
    invertAlgo.Det(d1, d2);
    *determinant = ldexp(d1, d2);
  }
}

} /* End of namespace genfit */
