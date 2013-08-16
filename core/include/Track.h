#ifndef genfit_Track_h
#define genfit_Track_h

#include <vector>

#include "AbsTrackRep.h"
#include "FitStatus.h"
#include "TrackCand.h"
#include "TrackPoint.h"

#include "TObject.h"
#include "TVectorD.h"

namespace genfit {

/**
 * Class for hit sorting.
 */
class TrackPointComparator {
 public:
  /**
   * Comparison operator used in Track::sort(). Compares sorting parameter.
   */
  bool operator() (const TrackPoint* lhs, const TrackPoint* rhs) const {
    return lhs->getSortingParameter() < rhs->getSortingParameter();
  }
};

  /** 
   *  Holds a number of #AbsTrackRep, which correspond to the different particle hypotheses or track models which should be fitted.
   *  
   *  Can be created from a TrackCand, then all the information from the TrackCand will go into the Track.
   */
class Track : public TObject {

 public:

  Track();
  Track(const TrackCand& trackCand);
  Track(AbsTrackRep* trackRep, const TVectorD& stateSeed);

  Track(const Track&); // copy constructor

  virtual ~Track();
  virtual void Clear(Option_t* = "");

  TrackPoint* getPoint(int id) const;
  const std::vector< genfit::TrackPoint* > & getPoints() const {return trackPoints_;}
  unsigned int getNumPoints() const {return trackPoints_.size();}

  TrackPoint* getPointWithMeasurement(int id) const;
  const std::vector< genfit::TrackPoint* > & getPointsWithMeasurement() const  {return trackPointsWithMeasurement_;}
  unsigned int getNumPointsWithMeasurement() const {return trackPointsWithMeasurement_.size();}

  TrackPoint* getPointWithMeasurementAndFitterInfo(int id, const AbsTrackRep* rep) const;

  AbsTrackRep* getTrackRep(int id) const {return trackReps_.at(id);}
  unsigned int getNumReps() const {return trackReps_.size();}

  // This is used when streaming TrackPoints.
  int getIdForRep(const AbsTrackRep* rep) const;

  /** @brief Get cardinal track representation
   *
   * The user has to choose which track rep should be considered the
   * best one after the fit. E.g. the track representation giving the
   * smallest chi2 could be chosen. By default the first in the list is returned.
   */
  AbsTrackRep* getCardinalRep() const {return trackReps_.at(cardinalRep_);}
  unsigned int getCardinalRepID() const {return cardinalRep_;}

  bool hasFitStatus(const AbsTrackRep* rep) const {if (fitStatuses_.find(rep) == fitStatuses_.end()) return false; return (fitStatuses_.at(rep) != NULL);}
  FitStatus* getFitStatus(const AbsTrackRep* rep) const {return fitStatuses_.at(rep);}
  void setFitStatus(FitStatus* fitStatus, const AbsTrackRep* rep);

  const TVectorD& getStateSeed() const {return stateSeed_;}
  void setStateSeed(const TVectorD& s) {stateSeed_.ResizeTo(s); stateSeed_ = s;}

  /** Insert TrackPoint BEFORE TrackPoint with position id, if id >= 0.
   * Id -1 means after last TrackPoint. Id -2 means before last TrackPoint. ...
   * Also deletes backwardInfos before new point and forwardInfos after new point.
   * Also sets Track backpointer of point accordingly.
   */
  void insertPoint(TrackPoint* point, int id = -1);

  void insertPoints(std::vector<genfit::TrackPoint*> points, int id = -1);

  void deletePoint(int id);

  /**
   * Merge two tracks. The TrackPoints of other will be inserted
   * after id (per default, they will be appended at the end).
   * The other track will not be altered, the TrackPoints will be (deep) copied.
   * Only copies the TrackPoints, NOT the reps, fit statuses and seed state.
   */
  void mergeTrack(const Track* other, int id = -1);

  void addTrackRep(AbsTrackRep* trackRep);

  /** Delete a #TrackRep and all corresponding #FitterInfos in the #TrackPoints
   */
  void deleteTrackRep(int id);
  void setCardinalRep(int id);

  /** Sort #TrackPoints and according to their sorting parameters.
   * Returns if the order of the trackPoints has actually changed.
   */
  bool sort();

  void deleteForwardInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = NULL); // delete in range [startId, endId]. If rep == NULL, delete for ALL reps, otherwise only for rep.
  void deleteBackwardInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = NULL); // delete in range [startId, endId]. If rep == NULL, delete for ALL reps, otherwise only for rep.
  void deleteReferenceInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = NULL); // delete in range [startId, endId]. If rep == NULL, delete for ALL reps, otherwise only for rep.
  void deleteMeasurementInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = NULL); // delete in range [startId, endId]. If rep == NULL, delete for ALL reps, otherwise only for rep.
  void deleteFitterInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = NULL); // delete in range [startId, endId]. If rep == NULL, delete for ALL reps, otherwise only for rep.

  //! get TrackLength between to trackPoints
  double getTrackLen(AbsTrackRep* rep, int startId = 0, int endId = -1) const;
  //! get time of flight in ns between to trackPoints
  double getTOF(AbsTrackRep* rep, int startId = 0, int endId = -1) const;

  void prune(const Option_t* = "");

  void Print(const Option_t* = "") const;

  bool checkConsistency() const;

 private:

  void trackHasChanged();

  Track& operator=(const Track&); // assignment operator  // delete until properly implemented

  void fillPointsWithMeasurement();

  std::vector<AbsTrackRep*> trackReps_; // Ownership
  unsigned int cardinalRep_; // THE selected rep, default = 0;

  std::vector<TrackPoint*> trackPoints_; // Ownership
  std::vector<TrackPoint*> trackPointsWithMeasurement_; //! helper

  std::map< const AbsTrackRep*, FitStatus* > fitStatuses_; // Ownership over FitStatus*

  TVectorD stateSeed_; // 6D: position, momentum

  ClassDef(Track,1)

};

} /* End of namespace genfit */

#endif // genfit_Track_h
