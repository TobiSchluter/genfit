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

#ifndef genfit_Tools_h
#define genfit_Tools_h

#include <TMatrixDSym.h>

namespace genfit {
namespace tools {

  /** @brief Invert a matrix, throwing GFException when inversion fails.
   * Optional calculation of determinant.
   */
  void invertMatrix(const TMatrixDSym& mat, TMatrixDSym& inv, double* determinant = NULL);
  /** @brief Same, replacing its argument.
   */
  void invertMatrix(TMatrixDSym& mat, double* determinant = NULL);

} /* End of namespace tools */
} /* End of namespace genfit */
/** @} */

#endif // genfit_Tools_h
