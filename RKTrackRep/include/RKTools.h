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

#include <array>

#ifndef genfit_RKTools_h
#define genfit_RKTools_h

namespace genfit {

// Array Matrix typedefs. They are needed for SSE optimization:
// gcc can vectorize loops only if the array sizes are known.
typedef std::array<double, 1*3> M1x3;
typedef std::array<double, 1*4> M1x4;
typedef std::array<double, 1*6> M1x6;
typedef std::array<double, 1*7> M1x7;
typedef std::array<double, 5*5> M5x5;
typedef std::array<double, 6*6> M6x6;
typedef std::array<double, 7*7> M7x7;
typedef std::array<double, 8*7> M8x7;
typedef std::array<double, 6*5> M6x5;
typedef std::array<double, 7*5> M7x5;
typedef std::array<double, 5*6> M5x6;
typedef std::array<double, 5*7> M5x7;

namespace RKTools {

  void J_pMTxcov5xJ_pM(const M5x7& J_pM, const M5x5& cov5, M7x7& out7);
  void J_pMTxcov5xJ_pM(const M5x6& J_pM, const M5x5& cov5, M6x6& out6);

  void J_MpTxcov7xJ_Mp(const M7x5& J_Mp, const M7x7& cov7, M5x5& out5);
  void J_MpTxcov6xJ_Mp(const M6x5& J_Mp, const M6x6& cov6, M5x5& out5);

  void J_MMTxcov7xJ_MM(const M7x7& J_MM, M7x7& cov7);

  void J_MMxJ_MM(M7x7& J_MM, const M7x7& J_MM_old);

  void J_pMxJ_MMxJ_Mp(const M5x7& J_pM, const M7x7& J_MM, const M7x5& J_Mp, M5x5& J_pp, bool MMproj = true);

  void printDim(const double* mat, unsigned int dimX, unsigned int dimY);

}

} /* End of namespace genfit */

#endif // genfit_RKTools_h

/** @} */
