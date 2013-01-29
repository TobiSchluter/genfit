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

#ifndef GFSTATEONPLANE_H
#define GFSTATEONPLANE_H

#include <vector>
#include <list>
#include <iostream>

#include <TObject.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <Math/ProbFunc.h>

#include "GFDetPlane.h"
#include "GFSmartPointers.h"

class TVector3;
class GFAbsRecoHit;



class GFStateOnPlane : public TObject {

 public:

  GFStateOnPlane();
  GFStateOnPlane(TVectorD& state, TMatrixDSym& cov);


  //! returns dimension of state vector
  unsigned int getDim() const {return fState.GetNrows();}

  inline const TVectorD& getState() const {return fState;}
  inline const TMatrixDSym& getCov() const {return fCov;}

  double getStateElem(int i) const {return fState(i);}
  double getCovElem(int i, int j) const {return fCov(i,j);}

  TVector3 getPos();
  TVector3 getMom();
  void getPosMom(TVector3& pos,TVector3& mom);

  //! method which gets position, momentum and 6x6 covariance matrix
  void getPosMomCov(TVector3& pos, TVector3& mom, TMatrixDSym& cov);

  //std::tr1::shared_ptr<GFDetPlane> getReferencePlane() const {return fPlane;}


  inline void setCov(const TMatrixDSym& cov) {fCov = cov;}

  //! method which sets position, momentum and 6x6 covariance matrix
  virtual void setPosMomCov(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov);


  //void Print(const Option_t* = "") const;

 protected:
  //! The vector of track parameters
  TVectorD fState;

  //! The covariance matrix
  TMatrixDSym fCov;

  // plane where the track parameters are given
  SHARED_PTR(GFDetPlane) fPlane;
  //GFDetPlane* fPlane;


 public:
  ClassDef(GFStateOnPlane,1)

};

#endif

/** @} */
