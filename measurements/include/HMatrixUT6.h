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

#ifndef genfit_HMatrixUT6_h
#define genfit_HMatrixUT6_h

#include "AbsHMatrix.h"


namespace genfit {

/**
 * @brief AbsHMatrix implementation for one-dimensional
 * MeasurementOnPlane with time information and RKTrackRepTime
 * parameterization.
 *
 * This projects out u and time.
 * H = (0, 0, 0, 1, 0, 0)
 *     (0, 0, 0, 0, 0, 1)
 */

class HMatrixUT6 : public AbsHMatrix {

 public:

  HMatrixUT6() {;}

  const TMatrixD& getMatrix() const;

  TVectorD Hv(const TVectorD& v) const;

  TMatrixD MHt(const TMatrixDSym& M) const;
  TMatrixD MHt(const TMatrixD& M) const;

  void HMHt(TMatrixDSym& M) const;

  virtual HMatrixUT6* clone() const {return new HMatrixUT6(*this);}

  virtual bool isEqual(const AbsHMatrix& other) const {return (dynamic_cast<const HMatrixUT6*>(&other) != NULL);}

  virtual void Print(const Option_t* = "") const;

  ClassDef(HMatrixUT6,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HMatrixUT6_h
