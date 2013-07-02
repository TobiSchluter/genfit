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

#include <algorithm>
#include <iostream>
#include <map>

//#include <glog/logging.h>


namespace genfit {

Track::Track() :
  cardinalRep_(0), stateSeed_(6)
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
  : cardinalRep_(rhs.cardinalRep_), stateSeed_(rhs.stateSeed_)
{
  std::map<const AbsTrackRep*, AbsTrackRep*> oldRepNewRep;

  for (std::vector<AbsTrackRep*>::const_iterator it=rhs.trackReps_.begin(); it!=rhs.trackReps_.end(); ++it) {
    AbsTrackRep* newRep = (*it)->clone();
    addTrackRep(newRep);
    oldRepNewRep[(*it)] = newRep;
  }

  for (std::vector<TrackPoint*>::const_iterator it=rhs.trackPoints_.begin(); it!=rhs.trackPoints_.end(); ++it) {
    insertPoint(new TrackPoint(**it, oldRepNewRep));
  }

  // self test
  assert(checkConsistency());
}

/*Track& Track::operator=(const Track& rhs) {
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

  assert(checkConsistency());

  return *this;
}*/


Track::~Track() {
  ; // smart pointers take care
}


TrackPoint* Track::getPoint(int id) const {
  if (id < 0)
    id += trackPoints_.size();

  return trackPoints_.at(id);
}


std::vector<TrackPoint*> Track::getPoints() const {
  std::vector< TrackPoint* > retVal;
  retVal.reserve(trackPoints_.size());

  for (std::vector<TrackPoint*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    retVal.push_back(*it);
  }

  return retVal;
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

  return NULL;
}


std::vector<TrackPoint*> Track::getPointsWithMeasurement() const {
  std::vector<TrackPoint*> retVal;

  for (std::vector<TrackPoint*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    if ((*it)->hasRawMeasurements()) {
      retVal.push_back(*it);
    }
  }

  return retVal;
}


unsigned int Track::getNumPointsWithMeasurement() const {
  unsigned int retVal(0);

  for (std::vector<TrackPoint*>::const_iterator it = trackPoints_.begin(); it != trackPoints_.end(); ++it) {
    if ((*it)->hasRawMeasurements()) {
      ++retVal;
    }
  }

  return retVal;
}


void Track::insertPoint(TrackPoint* point, int id) {
  // TODO: test
  point->setTrack(this);

  if (id == -1 || trackPoints_.size() == 0) {
    trackPoints_.push_back(point);
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

  trackPoints_.insert(trackPoints_.begin() + id + 1, point);  // insert inserts BEFORE
}


void Track::deletePoint(int id) {
  // TODO: test
  if (id < 0)
    id += trackPoints_.size();

  if (trackPoints_.at(id) != NULL) {
    // delete fitter infos if deleted point has a measurement
    if (trackPoints_[id]->hasRawMeasurements()) {
      deleteForwardInfo(id, -1);
      deleteBackwardInfo(0, id-1);
    }
  }
  trackPoints_.erase (trackPoints_.begin()+id);
}


void Track::mergeTrack(int i, Track other) {
  // TODO: implement
}


void Track::addTrackRep(AbsTrackRep* trackRep) {
  trackReps_.push_back( trackRep );
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
  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin(); pointIt != trackPoints_.end(); ++pointIt) {
    (*pointIt)->deleteFitterInfos(rep);
  }

  trackReps_.erase(trackReps_.begin()+id);
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

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    const std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
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

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    const std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
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

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
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

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
      (*fitterInfoIt)->deleteMeasurementInfo();
    }
  }
}


void Track::Print(const Option_t*) const {
  std::cout << "=======================================================================================\n";
  std::cout << "genfit::Track, containing " << trackPoints_.size() << " TrackPoints and " << trackReps_.size() << " TrackReps.\n";
  std::cout << " Seed state: "; stateSeed_.Print();

  for (unsigned int i=0; i<trackReps_.size(); ++i) {
    std::cout << " TrackRep Nr. " << i;
    if (i == cardinalRep_)
      std::cout << " (This is the cardinal rep)";
    std::cout << "\n";
    trackReps_[i]->Print();
  }

  std::cout << "---------------------------------------------------------------------------------------\n";

  for (unsigned int i=0; i<trackPoints_.size(); ++i) {
    std::cout << "TrackPoint Nr. " << i << "\n";
    trackPoints_[i]->Print();
    std::cout << "..........................................................................\n";
  }

  std::cout << "=======================================================================================\n";

}


bool Track::checkConsistency() const {
   
  // check if seed is 6D
  if (stateSeed_.GetNrows() != 6) {
    std::cerr << "Track::checkConsistency(): stateSeed_ dimension != 6" << std::endl;
    return false;
  }

  // check if cardinalRep_ is in range of trackReps_
  if (cardinalRep_ >= trackReps_.size()) {
    std::cerr << "Track::checkConsistency(): cardinalRep id out of bounds" << std::endl;
    return false;
  }

  // check TrackPoints
  for (std::vector<TrackPoint*>::const_iterator tp = trackPoints_.begin(); tp != trackPoints_.end(); ++tp) {
    // check if trackPoint points bach to this track
    if ((*tp)->getTrack() != this) {
      std::cerr << "Track::checkConsistency(): TrackPoint does not point back to this track" << std::endl;
      return false;
    }

    // check fitterInfos
    for (std::vector<AbsFitterInfo*>::const_iterator fi = ((*tp)->getFitterInfos()).begin(), fiend = ((*tp)->getFitterInfos()).end(); 
			fi < fiend; ++fi) {
      if (!( (*fi)->checkConsistency() ) ) {
        std::cerr << "Track::checkConsistency(): FitterInfo not consistent" << std::endl;
        return false;
      }

      // check if fitterInfos point to valid TrackReps in trackReps_
      int mycount (0);
      for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
        if ( (*rep) == (*fi)->getRep() ) {
          ++mycount;
        }
      }
      if (mycount ==  0) {
        std::cerr << "Track::checkConsistency(): fitterInfo points to TrackRep which is not in Track" << std::endl;
        return false;
      }
    }
  }

  return true;
}


} /* End of namespace genfit */
