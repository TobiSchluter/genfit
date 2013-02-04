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

#ifndef GFABSSTATEONPLANE_H
#define GFABSSTATEONPLANE_H

#include <TObject.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

#include "GFDetPlane.h"
#include "GFSmartPointers.h"


class GFAbsStateOnPlane : public TObject {

 public:

  GFAbsStateOnPlane(unsigned int dim = 0);
  GFAbsStateOnPlane(SHARED_PTR(GFDetPlane) plane, const TVectorD& state);

  //! return dimension of state
  unsigned int getDim() const {return fState.GetNrows();}

  const TVectorD& getState() const {return fState;}

  double getStateElem(int i) const {return fState(i);}

  const SHARED_PTR(GFDetPlane) getReferencePlane() const {return fPlane;}

  // accessors in cartesian coordinates
  virtual TVector3 getPos() const = 0;
  virtual TVector3 getMom() const = 0;
  virtual void getPosMom(TVector3& pos, TVector3& mom) const = 0;


  void setData(SHARED_PTR(GFDetPlane) plane, const TVectorD& state);
  void setState(const TVectorD& state);

  virtual void Print(const Option_t* = "") const;


 protected:
  //! plane where the track parameters are given
  SHARED_PTR(GFDetPlane) fPlane; // here a shared pointer is used, since different #GFAbsStateOnPlane objects can be defined in the same plane. Thus they can share ownership of that plane.

  //! The vector of track parameters
  TVectorD fState;


 public:
  ClassDef(GFAbsStateOnPlane,1)

};

#endif

/** @} */
