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

#include "SharedPlanePtr.h"


namespace genfit {

class AbsTrackRep;

  /** 
   *  A state with arbitrary dimension defined in a #GFDetPlane. #fSharedPlane is a shared_pointer, the ownership over that plane is shared between all #GFStateOnPlane objects defined in that plane.
   *  
   *  The definition of the state is bound to the TrackRep. Therefore, the #GFStateOnPlane contains a pointer to a #GFAbsTrackRep. It will provide functionality to extrapolate it and translate the state it into cartesian coordinates. 
   */
class StateOnPlane : public TObject {

 public:

  StateOnPlane(const AbsTrackRep* rep = _GFNULLPTR);
  StateOnPlane(const TVectorD& state, SharedPlanePtr plane, const AbsTrackRep* rep);


  const TVectorD& getState() const {return state_;}
  TVectorD& getState() {return state_;}
  const TVectorD& getAuxInfo() const {return auxInfo_;}
  TVectorD& getAuxInfo() {return auxInfo_;}
  SharedPlanePtr getPlane() const {return sharedPlane_;}
  const AbsTrackRep* getRep() const {return rep_;}

  void setState(const TVectorD& state) {if(state_.GetNrows() == 0) state_.ResizeTo(state); state_ = state;}
  void setPlane(SharedPlanePtr plane) {sharedPlane_ = plane;}
  void setStatePlane(const TVectorD& state, SharedPlanePtr plane) {state_ = state; sharedPlane_ = plane;}
  void setAuxInfo(const TVectorD& auxInfo) {if(auxInfo_.GetNrows() == 0) auxInfo_.ResizeTo(auxInfo); auxInfo_ = auxInfo;}
  void setRep(const AbsTrackRep* rep) {rep_ = rep;}

  virtual void Print(Option_t* option = "") const _GFOVERRIDE;

 protected:

  TVectorD state_; // state vector
  TVectorD auxInfo_; // auxiliary information (e.g. charge, flight direction etc.)
  SharedPlanePtr sharedPlane_; // Shared ownership.

 private:

  /** Pointer to TrackRep with respect to which StateOnPlane is defined
   */
  const AbsTrackRep* rep_; // No ownership


  //ClassDef(StateOnPlane,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_StateOnPlane_h
