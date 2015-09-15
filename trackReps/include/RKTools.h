/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_RKTools_h
#define genfit_RKTools_h

#include <algorithm>

namespace genfit {

/**
 * @brief a matrix class as a thin overlay over double[] providing
 * type-safety and some of the simplest operations for convenience.
 */
template <size_t nRows, size_t nCols>
struct RKMatrix {
  double vals[nRows * nCols];

  /**
   * @brief Row/column aware accessor.
   */
  double& operator()(size_t iRow, size_t iCol) {
    return vals[nCols*iRow + iCol];
  }
  /**
   * @brief Row/column read-only accessor.
   */
  const double& operator()(size_t iRow, size_t iCol) const {
    return vals[nCols*iRow + iCol];
  }
  /**
   * @brief Accessor to underlying array.
   */
  double& operator[](size_t n) {
    return vals[n];
  }
  /**
   * @brief Read-only accessor to underlying array.
   */
  const double& operator[](size_t n) const {
    return vals[n];
  }
  /**
   * @brief Iterator pointing towards beginning of underlying array.
   */
  double* begin() { return vals; }
  /**
   * @brief Iterator pointing past the last element of the underlying array.
   */
  double* end() { return vals + nRows * nCols; }
  /**
   * @brief Read-only iterator pointing towards beginning of underlying array.
   */
  const double* begin() const { return vals; }
  /**
   * @brief Read-only iterator pointing past the last element of the underlying array.
   */
  const double* end() const { return vals + nRows * nCols; }
  /**
   * @brief Matrix assignment operator.
   */
  RKMatrix<nRows, nCols>& operator=(const RKMatrix<nRows, nCols>& o) {
    std::copy(o.begin(), o.end(), this->begin());
    return *this;
  }
  /**
   * @brief Basic matrix addition operator +=.
   */
  RKMatrix<nRows, nCols>& operator+=(const RKMatrix<nRows, nCols>& o) {
    for (size_t i = 0; i < nRows; ++i)
      for (size_t j = 0; j < nCols; ++j)
        this->operator()(i, j) += o(i, j);
    return *this;
  }
  /**
   * @brief Write matrix to console.
   */
  void print() const;
};

/**
 * @brief Addition of two matrices.
 */
template<size_t nRows, size_t nCols>
RKMatrix<nRows, nCols> operator+(const RKMatrix<nRows, nCols>& left,
                                 const RKMatrix<nRows, nCols>& right)
{
  return (RKMatrix<nRows, nCols>(left) += right);
}

/**
 * @brief Addition of a scalar to a matrix.
 */
template<size_t nRows, size_t nCols>
RKMatrix<nRows, nCols> operator*(const double& left, const RKMatrix<nRows, nCols>& right)
{
  RKMatrix<nRows, nCols> result(right);
  for (size_t i = 0; i < nRows*nCols; ++i)
    result[i] *= left;
  return result;
}

//@{
/**
 * @brief Abbreviations for the various matrices used in the code.
 */
typedef RKMatrix<1, 3> M1x3;
typedef RKMatrix<1, 4> M1x4;
typedef RKMatrix<1, 7> M1x7;
typedef RKMatrix<1, 8> M1x8;
typedef RKMatrix<5, 5> M5x5;
typedef RKMatrix<6, 6> M6x6;
typedef RKMatrix<7, 7> M7x7;
typedef RKMatrix<8, 8> M8x8;
typedef RKMatrix<6, 5> M6x5;
typedef RKMatrix<7, 5> M7x5;
typedef RKMatrix<8, 5> M8x5;
typedef RKMatrix<8, 6> M8x6;
typedef RKMatrix<5, 6> M5x6;
typedef RKMatrix<5, 7> M5x7;
typedef RKMatrix<5, 8> M5x8;
typedef RKMatrix<6, 8> M6x8;
//@}

/**
 * @brief Array matrix multiplications used in RKTrackRep
 */
namespace RKTools {

  void J_pMTxcov5xJ_pM(const M5x7& J_pM, const M5x5& cov5, M7x7& out7);
  void J_pMTxcov5xJ_pM(const M5x6& J_pM, const M5x5& cov5, M6x6& out6);

  void J_pMTxcov6xJ_pM(const M6x8& J_pM, const M6x6& cov5, M8x8& out7);

  void J_MpTxnoise7xJ_Mp(const M8x6& J_Mp, const M7x7& noise7, M6x6& out6);
  void J_MpTxcov8xJ_Mp(const M8x6& J_Mp, const M8x8& cov8, M6x6& out6);
  void J_MpTxcov7xJ_Mp(const M7x5& J_Mp, const M7x7& cov7, M5x5& out5);
  void J_MpTxcov6xJ_Mp(const M6x5& J_Mp, const M6x6& cov6, M5x5& out5);

  void J_MMTxcov7xJ_MM(const M7x7& J_MM, M7x7& cov7);

  void J_MMxJ_MM(M7x7& J_MM, const M7x7& J_MM_old);

  void J_pMTTxJ_MMTTxJ_MpTT(const M8x6& J_pMT, const M8x8& J_MMT, const M6x8& J_MpT, M6x6& J_pp);
  void J_pMTTxJ_MMTTxJ_MpTT(const M7x5& J_pMT, const M7x7& J_MMT, const M5x7& J_MpT, M5x5& J_pp);

  void Np_N_NpT(const M7x7& Np, M7x7& N);

  void printDim(const double* mat, unsigned int dimX, unsigned int dimY);

}

/** @brief Implementation of the print function. */
template<size_t nRows, size_t nCols>
inline void
RKMatrix<nRows, nCols>::print() const {
  RKTools::printDim(this->vals, nRows, nCols);
}

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTools_h

