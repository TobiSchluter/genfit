#ifndef genfit_Track_h
#define genfit_Track_h

#include <vector>

#include "AbsTrackRep.h"
#include "Track.h"
#include "TrackCand.h"
#include "TrackPoint.h"


namespace genfit {


  /** 
   *  Holds a number of #AbsTrackRep, which correspond to the different particle hypotheses or track models which should be fitted.
   *  
   *  Can be created from a TrackCand, then all the information from the TrackCand will go into the Track.
   */
class Track {

 public:

  Track();
  Track(const TrackCand& trackCand);
  Track(AbsTrackRep* trackRep);

  Track(const Track&); // copy constructor
  Track(Track&&) = default; // move constructor
  Track& operator=(const Track&); // assignment operator
  Track& operator=(Track&&) = default; // move assignment operator

  ~Track();


  TrackPoint* getPoint(int id) const;
  const std::vector< TrackPoint* >& getPoints() {return trackPoints_;}
  unsigned int getNumPoints() const {return trackPoints_.size();}

  TrackPoint* getPointWithMeasurement(int id) const;
  std::vector< TrackPoint* > getPointsWithMeasurement() const;
  unsigned int getNumPointsWithMeasurement() const;

  AbsTrackRep* getTrackRep(int id) const {return trackReps_.at(id);}
  unsigned int getNumReps() const {return trackReps_.size();}

  /** @brief Get cardinal track representation
   *
   * The user has to choose which track rep should be considered the
   * best one after the fit. E.g. the track representation giving the
   * smallest chi2 could be chosen. By default the first in the list is returned.
   */
  AbsTrackRep* getCardinalRep() const {return trackReps_.at(cardinalRep_);}
  unsigned int getCardinalRepID() const {return cardinalRep_;}

  /** Insert TrackPoint before TrackPoint with position id.
   * Id -1 means after last TrackPoint.
   */
  void insertPoint(TrackPoint* point, int id = -1);

  void deletePoint(int id);

  void mergeTrack(int i, Track other);

  void addTrackRep(AbsTrackRep* trackRep);
  /** Delete a #TrackRep and all corresponding #FitterInfos in the #TrackPoints
   */
  void deleteTrackRep(int id);
  void setCardinalRep(int id);

  /** Sort #TrackPoints and according to their sorting parameters.
   */
  void sortHits();

  void deleteForwardInfo(int startId, int endId); // delete in range [startId, endId]
  void deleteBackwardInfo(int startId, int endId);
  void deleteReferenceInfo(int startId, int endId);
  void deleteMeasurementInfo(int startId, int endId);

 private:

  std::vector< TrackPoint* > trackPoints_; // Ownership

  std::vector< AbsTrackRep* > trackReps_; // Ownership
  unsigned int cardinalRep_; // THE selected rep, default = 0;


  //ClassDef(Track,1)

};

} /* End of namespace genfit */

#endif // genfit_Track_h
