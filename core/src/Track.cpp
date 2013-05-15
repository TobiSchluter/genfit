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

#include "Track.h"

#include <iostream>

//#include <glog/logging.h>


namespace genfit {

Track::Track() :
  cardinalRep_(0)
{
  ;
}

Track::Track(const TrackCand& trackCand) :
  cardinalRep_(0)
{
  // TODO: implement
}

Track::Track(AbsTrackRep* trackRep, const TVectorD& stateSeed) :
  cardinalRep_(0), stateSeed_(stateSeed)
{
  trackReps_.push_back(trackRep);
}


Track::Track(const Track& rhs)
  : stateSeed_(rhs.stateSeed_)
{
  for (TrackPoint* rhsPoint : rhs.trackPoints_) {
    trackPoints_.push_back(new TrackPoint(*rhsPoint));
  }

  for (AbsTrackRep* rhsRep : rhs.trackReps_) {
    trackReps_.push_back(rhsRep->clone());
  }

  cardinalRep_ = rhs.cardinalRep_;
}

Track& Track::operator=(const Track& rhs) {
  for (TrackPoint* point : trackPoints_) {
    if (point != nullptr)
      delete point;
  }
  trackPoints_.clear();

  for (TrackPoint* rhsPoint : rhs.trackPoints_) {
    trackPoints_.push_back(new TrackPoint(*rhsPoint));
  }

  for (AbsTrackRep* rep : trackReps_) {
    if (rep != nullptr)
      delete rep;
  }
  trackReps_.clear();

  for (AbsTrackRep* rhsRep : rhs.trackReps_) {
    trackReps_.push_back(rhsRep->clone());
  }

  cardinalRep_ = rhs.cardinalRep_;

  stateSeed_.ResizeTo(rhs.stateSeed_);
  stateSeed_ = rhs.stateSeed_;

  return *this;
}


Track::~Track() {
  for (TrackPoint* point : trackPoints_) {
    if (point != nullptr)
      delete point;
  }

  for (AbsTrackRep* rep : trackReps_) {
    if (rep != nullptr)
      delete rep;
  }
}


TrackPoint* Track::getPoint(int id) const {
  if (id < 0)
    id += trackPoints_.size();

  return trackPoints_.at(id);
}


TrackPoint* Track::getPointWithMeasurement(int id) const {
  if (id < 0)
    id += trackPoints_.size();

  int idMeas(0);

  for (unsigned int i=0; i<trackPoints_.size(); ++i) {
    if (trackPoints_[i]->hasRawMeasurements()) {
      if (id == idMeas)
        return trackPoints_.at(i);
      ++idMeas;
    }
  }

  return nullptr;
}


std::vector<TrackPoint*> Track::getPointsWithMeasurement() const {
  std::vector<TrackPoint*> retVal;

  for (unsigned int i=0; i<trackPoints_.size(); ++i) {
    if (trackPoints_[i]->hasRawMeasurements()) {
      retVal.push_back(trackPoints_[i]);
    }
  }

  return retVal;
}


unsigned int Track::getNumPointsWithMeasurement() const {
  unsigned int retVal(0);

  for (unsigned int i=0; i<trackPoints_.size(); ++i) {
    if (trackPoints_[i]->hasRawMeasurements()) {
      ++retVal;
    }
  }

  return retVal;
}


void Track::insertPoint(TrackPoint* point, int id) {
  // TODO: test

  if (trackPoints_.size() == 0) {
    assert(id == -1 || id == 0);
    trackPoints_.insert(trackPoints_.begin(), point);
    point->setTrack(this);
    return;
  }

  // [-size, size-1] is the allowed range
  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size();

  // delete fitter infos if inserted point has a measurement
  if (point->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id-1);
  }

  trackPoints_.insert(trackPoints_.begin() + id, point);
  point->setTrack(this);
}


void Track::deletePoint(int id) {
  // TODO: test
  if (id < 0)
    id += trackPoints_.size();

  if (trackPoints_.at(id) != nullptr) {
    // delete fitter infos if deleted point has a measurement
    if (trackPoints_[id]->hasRawMeasurements()) {
      deleteForwardInfo(id, -1);
      deleteBackwardInfo(0, id-1);
    }
    delete trackPoints_[id];
  }
  trackPoints_.erase (trackPoints_.begin()+id);
}


void Track::mergeTrack(int i, Track other) {
  // TODO: implement
}


void Track::addTrackRep(AbsTrackRep* trackRep) {
  // TODO: Test
  trackReps_.push_back(trackRep);

  // the fitter has to take care of adding fitterInfos
}


void Track::deleteTrackRep(int id) {
  // TODO: test
  if (id < 0)
    id += trackReps_.size();

  AbsTrackRep* rep = trackReps_.at(id);

  // update cardinalRep_
  if (int(cardinalRep_) == id)
    cardinalRep_ = 0; // reset
  else if (int(cardinalRep_) > id)
    --cardinalRep_; // make cardinalRep_ point to the same TrackRep before and after deletion

  // delete FitterInfos related to the deleted TrackRep
  for (auto pointIt = begin(trackPoints_); pointIt != end(trackPoints_); ++pointIt) {
    (*pointIt)->deleteFitterInfo(rep);
  }
}


void Track::setCardinalRep(int id) {
  // TODO: Test
  if (id < 0)
    id += trackReps_.size();

  if (id >= 0 && (unsigned int)id < trackReps_.size())
    cardinalRep_ = id;
  else {
    cardinalRep_ = 0;
    std::cerr << "Track::setCardinalRep: Attempted to set cardinalRep_ to a value out of bounds. Resetting  cardinalRep_ to 0." << std::endl;
  }
}


void Track::sortHits() {
  // TODO: implement
}


void Track::deleteForwardInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size() + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != end(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteForwardInfo();
    }
  }
}

void Track::deleteBackwardInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size() + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != end(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteBackwardInfo();
    }
  }
}

void Track::deleteReferenceInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size() + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != end(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteReferenceInfo();
    }
  }
}

void Track::deleteMeasurementInfo(int startId, int endId) {
  // TODO: test. Do we also have to delete forward- and backward info if measurements are removed?
  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size() + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != end(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteMeasurementInfo();
    }
  }
}


} /* End of namespace genfit */
