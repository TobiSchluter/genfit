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

#ifndef GFABSMEASUREDSTATEONPLANE_H
#define GFABSMEASUREDSTATEONPLANE_H

#include "GFAbsStateOnPlane.h"


class GFAbsMeasuredStateOnPlane : public GFAbsStateOnPlane {

 public:

  GFAbsMeasuredStateOnPlane(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, const TMatrixDSym& cov);
  GFAbsMeasuredStateOnPlane(const GFAbsStateOnPlane& state, const TMatrixDSym& cov);

  const TMatrixDSym& getCov() const {return fCov;}

  double getCovElem(int i, int j) const {return fCov(i,j);}


  //! method which gets position, momentum and 6x6 covariance matrix. Needed for #GFRave.
  virtual void getPosMomCov(TVector3& pos, TVector3& mom, TMatrixDSym& cov) = 0;


  void setData(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, const TMatrixDSym& cov);
  void setStateCov(const TVectorD& state, const TMatrixDSym& cov);
  void setCov(const TMatrixDSym& cov);

  //! method which sets position, momentum and 6x6 covariance matrix. Needed for #GFRave.
  virtual void setPosMomCov(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov) = 0;

  virtual void Print(const Option_t* = "") const;


 protected:

  void checkDim();

  //! The covariance matrix
  TMatrixDSym fCov;

 public:
  ClassDef(GFAbsMeasuredStateOnPlane,1)

};

#endif

/** @} */
