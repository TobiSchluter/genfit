/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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
  Track(AbsTrackRep* trackRep, const TVectorD& stateSeed, const TMatrixDSym& covSeed);

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

  /** Shortcut to get FittedStates. Uses getPointWithMeasurementAndFitterInfo(id, rep).
   * Per default, the fitted state of the fitterInfo of the first TrackPoint with measurements(s) and fitterInfo(s)
   * is returned. If no trackRep is specified, the fitter info of the cardinal rep will be used.
   */
  const MeasuredStateOnPlane& getFittedState(int id = 0, const AbsTrackRep* rep = NULL, bool biased = true) const;

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
  unsigned int getCardinalRepId() const {return cardinalRep_;}

  //! Check if track has a fitStatus for given rep. Per default, check for cardinal rep.
  bool hasFitStatus(const AbsTrackRep* rep = NULL) const;
  //! Get fit status for a TrackRep. Per default, return FitStatus for cardinalRep.
  FitStatus* getFitStatus(const AbsTrackRep* rep = NULL) const {if (rep == NULL) rep = getCardinalRep(); return fitStatuses_.at(rep);}
  void setFitStatus(FitStatus* fitStatus, const AbsTrackRep* rep);

  const TVectorD& getStateSeed() const {return stateSeed_;}
  void setStateSeed(const TVectorD& s) {stateSeed_.ResizeTo(s); stateSeed_ = s;}

  const TMatrixDSym& getCovSeed() const {return covSeed_;}
  void setCovSeed(const TMatrixDSym& c) {covSeed_.ResizeTo(c); covSeed_ = c;}

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

  /**
   * C:  prune all reps except cardinalRep
   * F:  prune all points except first point
   * L:  prune all points except last point
   * FL: prune all points except first and last point
   * W:  prune rawMeasurements from TrackPoints
   * R:  prune referenceInfo from fitterInfos
   * M:  prune measurementInfo from fitterInfos
   * I:  if F, L, or FL is set, prune forward (backward) info of first (last) point
   * U:  if fitterInfo is a KalmanFitterInfo, prune predictions and keep updates
   */
  void prune(const Option_t* = "CFLWRMIU");

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
  TMatrixDSym covSeed_; // 6D

  ClassDef(Track,1)

};

} /* End of namespace genfit */

#endif // genfit_Track_h
