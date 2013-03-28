#include "StateOnPlane.h"

namespace genfit {


  /** 
   *  A state with arbitrary dimension defined in a #GFDetPlane. #fSharedPlane is a shared_pointer, the ownership over that plane is shared between all #GFStateOnPlane objects defined in that plane.
   *  
   *  The definition of the state is bound to the TrackRep. Therefore, the #GFStateOnPlane contains a pointer to a #GFAbsTrackRep. It will provide functionality to extrapolate it and translate the state it into cartesian coordinates. 
   *  If the #GFStateOnPlane is calculated from a #GFAbsRawMeasurement, it does in general not have the full dimensionality and can therefore not be extrapolated in general.
   */



} /* End of namespace genfit */
