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

#ifndef genfit_MeasuredStateOnPlane_h
#define genfit_MeasuredStateOnPlane_h

#include <TMatrixDSym.h>

#include "StateOnPlane.h"


namespace genfit {


  /** 
   *  Additional covariance matrix.
   */
class MeasuredStateOnPlane : public StateOnPlane {

 public:

  MeasuredStateOnPlane(const AbsTrackRep* rep = nullptr);
  MeasuredStateOnPlane(const TVectorD& state, const TMatrixDSym& cov, SharedPlanePtr plane, const AbsTrackRep* rep);
  MeasuredStateOnPlane(const StateOnPlane& state, const TMatrixDSym& cov);

  const TMatrixDSym& getCov() const {return cov_;}
  TMatrixDSym& getCov() {return cov_;}

  void setStateCov(const TVectorD& state, const TMatrixDSym& cov) {setState(state); setCov(cov);}
  void setStateCovPlane(const TVectorD& state, const TMatrixDSym& cov, SharedPlanePtr plane) {setStatePlane(state, plane); setCov(cov);}
  void setCov(const TMatrixDSym& cov) {if(cov_.GetNrows() == 0) cov_.ResizeTo(cov); cov_ = cov;}

  virtual void Print(Option_t* option = "") const override;

 protected:

  TMatrixDSym cov_;


  //ClassDef(MeasuredStateOnPlane,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_MeasuredStateOnPlane_h
