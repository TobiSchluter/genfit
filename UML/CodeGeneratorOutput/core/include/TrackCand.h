#ifndef genfit_TrackCand_h
#define genfit_TrackCand_h

#include <vector>

#include "TrackCandHit.h"


namespace genfit {

class TrackCand {

 public:

  /**
   * @element-type TrackCandHit
   */
  std::vector< TrackCandHit > myTrackCandHit;
};

} /* End of namespace genfit */

#endif // genfit_TrackCand_h
