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


void Track::removePoint(unsigned int i) {
  // TODO: implement
}


void Track::mergeTrack(int i, Track other) {
  // TODO: implement
}


void Track::addTrackRep(AbsTrackRep* trackRep) {
  // TODO: implement
}


void Track::removeTrackRep(unsigned int i) {
  // TODO: implement
}


void Track::setCardinalRep(unsigned int id) {
  // TODO: implement
}


void Track::sortHits() {
  // TODO: implement
}


} /* End of namespace genfit */
