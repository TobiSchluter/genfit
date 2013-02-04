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

#ifndef GFSTATEONPLANE5D_H
#define GFSTATEONPLANE5D_H

#include "GFAbsStateOnPlane.h"

#define DIM 5


/** @brief The state is 5D: (q/p, u', v', u, v)
 */

class GFStateOnPlane5D : virtual public GFAbsStateOnPlane {

 public:

  GFStateOnPlane5D();
  GFStateOnPlane5D(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, double charge, bool spu);

  void setData(SHARED_PTR(GFDetPlane) plane, const TVectorD& state, double charge, bool spu);

  // accessors in cartesian coordinates
  TVector3 getPos() const;
  TVector3 getMom() const;
  void getPosMom(TVector3& pos, TVector3& mom) const;

  double getCharge() const {return fCharge;}
  bool getSpu() const {return fSpu;}

 protected:

  //! electric charge of the particle
  double fCharge;

  //! flight direction. true: Mom (dot) fPlane.getNormal() > 0; false: Mom (dot) fPlane.getNormal() < 0
  bool fSpu;


 public:
  ClassDef(GFStateOnPlane5D,1)

};

#endif

/** @} */
