#include "StateOnPlane.h"

namespace genfit {


  /** 
   *  A state with arbitrary dimension defined in a #GFDetPlane. #fSharedPlane is a shared_pointer, the ownership over that plane is shared between all #GFStateOnPlane objects defined in that plane.
   *  
   *  The definition of the state is bound to the TrackRep. Therefore, the #GFStateOnPlane contains a pointer to a #GFAbsTrackRep. It will provide functionality to extrapolate it and translate the state it into cartesian coordinates. 
   *  If the #GFStateOnPlane is calculated from a #GFAbsRawMeasurement, it does in general not have the full dimensionality and can therefore not be extrapolated in general.
   */



StateOnPlane::StateOnPlane(DetPlane plane, TVectorD state)
// don't delete the following line as it's needed to preserve source code of this autogenerated element
// section 127-0-1-1--1eeb8d28:13caa78033d:-8000:00000000000011B4 begin
{
}
// section 127-0-1-1--1eeb8d28:13caa78033d:-8000:00000000000011B4 end
// don't delete the previous line as it's needed to preserve source code of this autogenerated element

} /* End of namespace genfit */
