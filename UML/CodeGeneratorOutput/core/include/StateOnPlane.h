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

#ifndef genfit_StateOnPlane_h
#define genfit_StateOnPlane_h

#include <TObject.h>
#include <TVectorD.h>

#include "DetPlane.h"

namespace genfit {
class AbsTrackRep;
} /* End of namespace genfit */

namespace genfit {


  /** 
   *  A state with arbitrary dimension defined in a #GFDetPlane. #fSharedPlane is a shared_pointer, the ownership over that plane is shared between all #GFStateOnPlane objects defined in that plane.
   *  
   *  The definition of the state is bound to the TrackRep. Therefore, the #GFStateOnPlane contains a pointer to a #GFAbsTrackRep. It will provide functionality to extrapolate it and translate the state it into cartesian coordinates. 
   */
class StateOnPlane : public TObject {

 public:

  StateOnPlane();
  StateOnPlane(const TVectorD& state, DetPlane* plane, AbsTrackRep* rep);

  ~StateOnPlane();

  const TVectorD& getState() const {return state_;}
  const DetPlane* getPlane() const {return sharedPlane_;}
  const AbsTrackRep* getRep() const {return rep_;}

 protected:

  void setState(const TVectorD& state) {state_ = state;}
  void setStatePlane(const TVectorD& state, DetPlane* plane);

 protected:

  TVectorD state_;
  DetPlane* sharedPlane_; // Ownership. TODO: Change to std::shared_ptr when changing to ROOT 6
  AbsTrackRep* rep_; // No ownership


  ClassDef(StateOnPlane,1)

};

} /* End of namespace genfit */

#endif // genfit_StateOnPlane_h
