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

#ifndef genfit_KalmanFittedStateOnPlane_h
#define genfit_KalmanFittedStateOnPlane_h

#include "MeasuredStateOnPlane.h"


namespace genfit {


  /** 
   *  Additional info produced by a Kalman filter or DAF.
   */
class KalmanFittedStateOnPlane : public MeasuredStateOnPlane {

 public:

  KalmanFittedStateOnPlane();
  KalmanFittedStateOnPlane(const TVectorD& state, const TMatrixDSym& cov, SharedPlanePtr plane, const AbsTrackRep* rep, double chiSquareIncrement, double ndf);
  KalmanFittedStateOnPlane(const MeasuredStateOnPlane& state, double chiSquareIncrement, double ndf);

  double getChiSquareIncrement() const {return chiSquareIncrement_;}
  double getNdf() const {return ndf_;}

  void setChiSquareIncrement(double chiSquareIncrement) {chiSquareIncrement_ = chiSquareIncrement;}
  void setNdf(double ndf) {ndf_ = ndf;}


 protected:

  double chiSquareIncrement_;
  
  /** 
   *  Degrees of freedom. Needs to be a double because of DAF.
   */
  double ndf_;


  //ClassDef(KalmanFittedStateOnPlane,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFittedStateOnPlane_h
