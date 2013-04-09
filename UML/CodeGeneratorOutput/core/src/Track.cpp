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

Track::Track(AbsTrackRep* trackRep) :
  cardinalRep_(0)
{
  trackReps_.push_back(trackRep);
}


Track::Track(const Track& rhs) {
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


void Track::insertPoint(TrackPoint* point, int i) {
  // TODO: implement
}


void Track::deletePoint(int id) {
  if (id < 0)
    id = trackPoints_.size() + id;

  if (trackPoints_.at(id) != nullptr)
    delete trackPoints_[id];
  trackPoints_.erase (trackPoints_.begin()+id);


}


void Track::mergeTrack(int i, Track other) {
  // TODO: implement
}


void Track::addTrackRep(AbsTrackRep* trackRep) {
  // TODO: implement
}


void Track::deleteTrackRep(int id) {
  // TODO: implement
}


void Track::setCardinalRep(unsigned int id) {
  // TODO: implement
}


void Track::sortHits() {
  // TODO: implement
}


void Track::deleteForwardInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId = trackPoints_.size() + startId;
  if (endId < 0)
    endId = trackPoints_.size() + endId + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != begin(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteForwardInfo();
    }
  }
}

void Track::deleteBackwardInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId = trackPoints_.size() + startId;
  if (endId < 0)
    endId = trackPoints_.size() + endId + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != begin(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteBackwardInfo();
    }
  }
}

void Track::deleteReferenceInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId = trackPoints_.size() + startId;
  if (endId < 0)
    endId = trackPoints_.size() + endId + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != begin(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteReferenceInfo();
    }
  }
}

void Track::deleteMeasurementInfo(int startId, int endId) {
  // TODO: test
  if (startId < 0)
    startId = trackPoints_.size() + startId;
  if (endId < 0)
    endId = trackPoints_.size() + endId + 1;

  for (auto pointIt = begin(trackPoints_) + startId; pointIt != begin(trackPoints_) + endId; ++pointIt) {
    auto fitterInfos = (*pointIt)->getFitterInfos();
    for (auto fitterInfoIt = begin(fitterInfos); fitterInfoIt != begin(fitterInfos); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteMeasurementInfo();
    }
  }
}


} /* End of namespace genfit */
