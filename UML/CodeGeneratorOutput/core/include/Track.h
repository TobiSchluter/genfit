#ifndef genfit_Track_h
#define genfit_Track_h

#include <vector>

#include "AbsTrackRep.h"
#include "Track.h"
#include "TrackPoint.h"


namespace genfit {


  /** 
   *  Holds a number of #GFAbsTrackRep, which correspond to the different particle hypotheses or track models which should be fitted.
   *  
   *  Can be created from a TrackCand, then all the information from the TrackCand will go into the Track.
   */
class Track {

 public:

  virtual insertPoint(int i, TrackPoint point);

  virtual removePoint(int i);

  virtual mergeTrack(int i, Track other);

 public:

  /**
   * @element-type TrackPoint
   */
  std::vector< TrackPoint > trackPoints_;

  /**
   * @element-type AbsTrackRep
   */
  std::vector< AbsTrackRep* > trackReps_;
};

} /* End of namespace genfit */

#endif // genfit_Track_h
