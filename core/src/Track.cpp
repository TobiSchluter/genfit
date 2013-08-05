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
#include "KalmanFitterInfo.h"

#include <algorithm>
#include <iostream>
#include <map>

#include <TDatabasePDG.h>

//#include <glog/logging.h>

//#define DEBUG


namespace genfit {

Track::Track() :
  cardinalRep_(0), fitStatuses_(), stateSeed_(6)
{
  ;
}

Track::Track(const TrackCand& trackCand) :
  cardinalRep_(0)
{
  // TODO: implement
}

Track::Track(AbsTrackRep* trackRep, const TVectorD& stateSeed) :
  cardinalRep_(0), fitStatuses_(), stateSeed_(stateSeed)
{
  addTrackRep(trackRep);
}


Track::Track(const Track& rhs)
  : cardinalRep_(rhs.cardinalRep_), stateSeed_(rhs.stateSeed_)
{
  assert(rhs.checkConsistency());

  std::map<const AbsTrackRep*, AbsTrackRep*> oldRepNewRep;

  for (std::vector<AbsTrackRep*>::const_iterator it=rhs.trackReps_.begin(); it!=rhs.trackReps_.end(); ++it) {
    AbsTrackRep* newRep = (*it)->clone();
    addTrackRep(newRep);
    oldRepNewRep[(*it)] = newRep;
  }

  for (std::vector<TrackPoint*>::const_iterator it=rhs.trackPoints_.begin(); it!=rhs.trackPoints_.end(); ++it) {
    trackPoints_.push_back(new TrackPoint(**it, oldRepNewRep));
    trackPoints_.back()->setTrack(this);
  }

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=rhs.fitStatuses_.begin(); it!=rhs.fitStatuses_.end(); ++it) {
    setFitStatus(it->second->clone(), oldRepNewRep[it->first]);
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
  // FIXME: smarter containers or pointers needed ...
  for (size_t i = 0; i < trackPoints_.size(); ++i)
    delete trackPoints_[i];

  for (std::map< const AbsTrackRep*, FitStatus* >::iterator it = fitStatuses_.begin(); it!= fitStatuses_.end(); ++it)
    delete it->second;

  for (size_t i = 0; i < trackReps_.size(); ++i)
    delete trackReps_[i];
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


int Track::getIdForRep(const AbsTrackRep* rep) const
{
  for (size_t i = 0; i < trackReps_.size(); ++i)
    if (trackReps_[i] == rep)
      return i;

  assert(0 == 1);  // Cannot happen.
}


void Track::setFitStatus(FitStatus* fitStatus, const AbsTrackRep* rep) {
  if (fitStatuses_.find(rep) != fitStatuses_.end())
    delete fitStatuses_.at(rep);

  fitStatuses_[rep] = fitStatus;
}


void Track::insertPoint(TrackPoint* point, int id) {

  point->setTrack(this);

  #ifdef DEBUG
  std::cout << "Track::insertPoint at position " << id  << "\n";
  #endif
  assert(point!=NULL);
  trackHasChanged();

  point->setTrack(this);

  if (trackPoints_.size() == 0) {
    trackPoints_.push_back(point);
    return;
  }

  if (id == -1 || id == (int)trackPoints_.size()) {
    trackPoints_.push_back(point);
    deleteReferenceInfo(std::max(0, (int)trackPoints_.size()-2), (int)trackPoints_.size()-1);

    // delete fitter infos if inserted point has a measurement
    if (point->hasRawMeasurements()) {
      deleteForwardInfo(-1, -1);
      deleteBackwardInfo(0, -2);
    }

    return;
  }

  // [-size, size-1] is the allowed range
  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size() + 1;

  // insert
  trackPoints_.insert(trackPoints_.begin() + id, point);  // insert inserts BEFORE

  // delete fitter infos if inserted point has a measurement
  if (point->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id);
  }

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id+1));
}


void Track::insertPoints(std::vector<genfit::TrackPoint*> points, int id) {

  int nBefore = getNumPoints();
  int n = points.size();

  if (n == 0)
    return;
  if (n == 1)
    insertPoint(points[0], id);

  for (std::vector<genfit::TrackPoint*>::iterator p = points.begin(); p != points.end(); ++p)
    (*p)->setTrack(this);

  if (id == -1 || id == (int)trackPoints_.size()) {
    trackPoints_.insert(trackPoints_.end(), points.begin(), points.end());

    deleteReferenceInfo(std::max(0, nBefore-1), nBefore);

    deleteForwardInfo(nBefore, -1);
    deleteBackwardInfo(0, std::max(0, nBefore-1));

    return;
  }


  assert(id < (ssize_t)trackPoints_.size() || -id-1 <= (ssize_t)trackPoints_.size());

  if (id < 0)
    id += trackPoints_.size() + 1;


  // insert
  trackPoints_.insert(trackPoints_.begin() + id, points.begin(), points.end());  // insert inserts BEFORE

  // delete fitter infos if inserted point has a measurement
  deleteForwardInfo(id, -1);
  deleteBackwardInfo(0, id+n);

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id));
  deleteReferenceInfo(std::max(0, id+n-1), std::min((int)trackPoints_.size()-1, id+n));

}


void Track::deletePoint(int id) {

  #ifdef DEBUG
  std::cout << "Track::deletePoint at position " << id  << "\n";
  #endif

  // TODO: test
  trackHasChanged();

  if (id < 0)
    id += trackPoints_.size();
  assert(id>0);


  // delete forwardInfo after point (backwardInfo before point) if deleted point has a measurement
  if (trackPoints_[id]->hasRawMeasurements()) {
    deleteForwardInfo(id, -1);
    deleteBackwardInfo(0, id-1);
  }

  // delete reference info of neighbouring points
  deleteReferenceInfo(std::max(0, id-1), std::min((int)trackPoints_.size()-1, id+1));

  // delete point
  delete trackPoints_[id];
  trackPoints_.erase(trackPoints_.begin()+id);

}


void Track::mergeTrack(const Track* other, int id) {

  if (other->getNumPoints() == 0)
    return;

  std::map<const AbsTrackRep*, AbsTrackRep*> thisRepOtherRep;
  std::vector<const AbsTrackRep*> otherRepsToRemove;

  for (std::vector<AbsTrackRep*>::const_iterator thisRep=trackReps_.begin(); thisRep!=trackReps_.end(); ++thisRep) {
    for (std::vector<AbsTrackRep*>::const_iterator otherRep=other->trackReps_.begin(); otherRep!=other->trackReps_.end(); ++otherRep) {
      if ((*thisRep)->isSame(*otherRep))
        thisRepOtherRep[(*thisRep)] = *otherRep;
      else
        otherRepsToRemove.push_back(*otherRep);
    }
  }

  std::vector<TrackPoint*> points;
  points.reserve(other->getNumPoints());

  for (std::vector<TrackPoint*>::const_iterator otherTp=other->trackPoints_.begin(); otherTp!=other->trackPoints_.end(); ++otherTp) {
    points.push_back(new TrackPoint(**otherTp, thisRepOtherRep, &otherRepsToRemove));
  }

  insertPoints(points, id);
}


void Track::addTrackRep(AbsTrackRep* trackRep) {
  trackReps_.push_back(trackRep);
  fitStatuses_[trackRep] = new FitStatus();
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
    (*pointIt)->deleteFitterInfo(rep);
  }

  // delete fitStatus
  delete fitStatuses_.at(rep);
  fitStatuses_.erase(rep);

  // delete rep
  delete rep;
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


bool Track::sort() {
  #ifdef DEBUG
  std::cout << "Track::sort \n";
  #endif

  int nPoints(trackPoints_.size());
  // original order
  std::vector<TrackPoint*> pointsBefore(trackPoints_);

  // sort
  std::stable_sort(trackPoints_.begin(), trackPoints_.end(), TrackPointComparator());

  // see where order changed
  int equalUntil(-1), equalFrom(nPoints);
  for (int i = 0; i<nPoints; ++i) {
    if (pointsBefore[i] == trackPoints_[i])
      equalUntil = i;
    else
      break;
  }

  if (equalUntil == nPoints-1)
    return false; // sorting did not change anything


  trackHasChanged();

  for (int i = nPoints-1; i>equalUntil; --i) {
    if (pointsBefore[i] == trackPoints_[i])
      equalFrom = i;
    else
      break;
  }

  #ifdef DEBUG
  std::cout << "Track::sort. Equal up to (including) hit " << equalUntil << " and from (including) hit " << equalFrom << " \n";
  #endif

  deleteForwardInfo(equalUntil+1, -1);
  deleteBackwardInfo(0, equalFrom-1);

  deleteReferenceInfo(std::max(0, equalUntil+1), std::min((int)trackPoints_.size()-1, equalFrom-1));

  return true;
}


void Track::deleteForwardInfo(int startId, int endId, const AbsTrackRep* rep) {
  #ifdef DEBUG
  std::cout << "Track::deleteForwardInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != NULL) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteForwardInfo();
    }
    else {
      const std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteForwardInfo();
      }
    }
  }
}

void Track::deleteBackwardInfo(int startId, int endId, const AbsTrackRep* rep) {
  // TODO: test

  #ifdef DEBUG
  std::cout << "Track::deleteBackwardInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);


  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != NULL) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteBackwardInfo();
    }
    else {
      const std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteBackwardInfo();
      }
    }
  }
}

void Track::deleteReferenceInfo(int startId, int endId, const AbsTrackRep* rep) {
  // TODO: test

  #ifdef DEBUG
  std::cout << "Track::deleteReferenceInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != NULL) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteReferenceInfo();
    }
    else {
      std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteReferenceInfo();
      }
    }
  }
}

void Track::deleteMeasurementInfo(int startId, int endId, const AbsTrackRep* rep) {

  // TODO: test. Do we also have to delete forward- and backward info if measurements are removed?
  #ifdef DEBUG
  std::cout << "Track::deleteMeasurementInfo from position " << startId  << " to " << endId << "\n";
  #endif

  trackHasChanged();

  if (startId < 0)
    startId += trackPoints_.size();
  if (endId < 0)
    endId += trackPoints_.size();
  endId += 1;

  assert (endId >= startId);

  for (std::vector<TrackPoint*>::const_iterator pointIt = trackPoints_.begin() + startId; pointIt != trackPoints_.begin() + endId; ++pointIt) {
    if (rep != NULL) {
      if ((*pointIt)->hasFitterInfo(rep))
        (*pointIt)->getFitterInfo(rep)->deleteMeasurementInfo();
    }
    else {
      std::vector<AbsFitterInfo*> fitterInfos = (*pointIt)->getFitterInfos();
      for (std::vector<AbsFitterInfo*>::const_iterator fitterInfoIt = fitterInfos.begin(); fitterInfoIt != fitterInfos.end(); ++fitterInfoIt) {
        (*fitterInfoIt)->deleteMeasurementInfo();
      }
    }
  }
}


void Track::Print(const Option_t* option) const {
  TString opt = option;
  opt.ToUpper();
  if (opt.Contains("C")) { // compact

    std::cout << "\n   ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {

      int color = 32*(size_t)(trackPoints_[i]) % 15;
      switch (color) {
        case 0:
          std::cout<<"\033[1;30m";
          break;
        case 1:
          std::cout<<"\033[0;34m";
          break;
        case 2:
          std::cout<<"\033[1;34m";
          break;
        case 3:
          std::cout<<"\033[0;32m";
          break;
        case 4:
          std::cout<<"\033[1;32m";
          break;
        case 5:
          std::cout<<"\033[0;36m";
          break;
        case 6:
          std::cout<<"\033[1;36m";
          break;
        case 7:
          std::cout<<"\033[0;31m";
          break;
        case 8:
          std::cout<<"\033[1;31m";
          break;
        case 9:
          std::cout<<"\033[0;35m";
          break;
        case 10:
          std::cout<<"\033[1;35m";
          break;
        case 11:
          std::cout<<"\033[0;33m";
          break;
        case 12:
          std::cout<<"\033[1;33m";
          break;
        case 13:
          std::cout<<"\033[0;37m";
          break;
        default:
          ;
      }
      std::cout << trackPoints_[i] << "\033[00m  ";
    }
    std::cout << "\n";

    std::cout << "  ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {
      printf("% -9.3g  ", trackPoints_[i]->getSortingParameter());
    }
    std::cout << "\n";

    std::cout << "   ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {
      if (! trackPoints_[i]->hasFitterInfo(getCardinalRep())) {
        std::cout << "           ";
        continue;
      }
      AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(getCardinalRep());
      if (fi->hasMeasurements())
        std::cout << "M";
      else
        std::cout << " ";

      if (fi->hasReferenceState())
        std::cout << "R";
      else
        std::cout << " ";

      std::cout << "         ";
    }
    std::cout << "\n";

    std::cout << "-> ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {
      if (! trackPoints_[i]->hasFitterInfo(getCardinalRep())) {
        std::cout << "           ";
        continue;
      }
      AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(getCardinalRep());
      if (fi->hasForwardPrediction())
        std::cout << "P";
      else
        std::cout << " ";

      if (fi->hasForwardUpdate())
        std::cout << "U";
      else
        std::cout << " ";

      std::cout << "         ";
    }
    std::cout << "\n";

    std::cout << "<- ";
    for (unsigned int i=0; i<trackPoints_.size(); ++i) {
      if (! trackPoints_[i]->hasFitterInfo(getCardinalRep())) {
        std::cout << "           ";
        continue;
      }
      AbsFitterInfo* fi = trackPoints_[i]->getFitterInfo(getCardinalRep());
      if (fi->hasBackwardPrediction())
        std::cout << "P";
      else
        std::cout << " ";

      if (fi->hasBackwardUpdate())
        std::cout << "U";
      else
        std::cout << " ";

      std::cout << "         ";
    }
    std::cout << "\n";

    return;
  }



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

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->Print();
  }

  std::cout << "=======================================================================================\n";

}


bool Track::checkConsistency() const {

  std::map<const AbsTrackRep*, double> prevSegmentLengths;

  // check if seed is 6D
  if (stateSeed_.GetNrows() != 6) {
    std::cerr << "Track::checkConsistency(): stateSeed_ dimension != 6" << std::endl;
    return false;
  }

  // check if cardinalRep_ is in range of trackReps_
  if (trackReps_.size() && cardinalRep_ >= trackReps_.size()) {
    std::cerr << "Track::checkConsistency(): cardinalRep id " << cardinalRep_ << " out of bounds" << std::endl;
    return false;
  }

  for (std::vector<AbsTrackRep*>::const_iterator rep = trackReps_.begin(); rep != trackReps_.end(); ++rep) {
    // check for NULL
    if ((*rep) == NULL) {
      std::cerr << "Track::checkConsistency(): TrackRep is NULL" << std::endl;
      return false;
    }

    // check for valid pdg code
    TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle((*rep)->getPDG());
    if (particle == NULL) {
      std::cerr << "Track::checkConsistency(): TrackRep pdg ID is not valid" << std::endl;
      return false;
    }

    prevSegmentLengths[*rep] = 0;
  }

  // check TrackPoints
  for (std::vector<TrackPoint*>::const_iterator tp = trackPoints_.begin(); tp != trackPoints_.end(); ++tp) {
    // check for NULL
    if ((*tp) == NULL) {
      std::cerr << "Track::checkConsistency(): TrackPoint is NULL" << std::endl;
      return false;
    }
    // check if trackPoint points back to this track
    if ((*tp)->getTrack() != this) {
      std::cerr << "Track::checkConsistency(): TrackPoint does not point back to this track" << std::endl;
      return false;
    }

    // check rawMeasurements
    std::vector<AbsMeasurement*> rawMeasurements = (*tp)->getRawMeasurements();
    for (std::vector<AbsMeasurement*>::const_iterator m = rawMeasurements.begin(); m != rawMeasurements.end(); ++m) {
      // check for NULL
      if ((*m) == NULL) {
        std::cerr << "Track::checkConsistency(): Measurement is NULL" << std::endl;
        return false;
      }
      // check if measurement points back to TrackPoint
      if ((*m)->getTrackPoint() != *tp) {
        std::cerr << "Track::checkConsistency(): Measurement does not point back to correct TrackPoint" << std::endl;
        return false;
      }
    }

    // check fitterInfos
    std::vector<AbsFitterInfo*> fitterInfos = (*tp)->getFitterInfos();
    for (std::vector<AbsFitterInfo*>::const_iterator fi = fitterInfos.begin(); fi != fitterInfos.end(); ++fi) {
      // check for NULL
      if ((*fi) == NULL) {
        std::cerr << "Track::checkConsistency(): FitterInfo is NULL. TrackPoint: " << *tp << std::endl;
        return false;
      }

      if (!( (*fi)->checkConsistency() ) ) {
        std::cerr << "Track::checkConsistency(): FitterInfo not consistent. TrackPoint: " << *tp << std::endl;
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

      if (dynamic_cast<KalmanFitterInfo*>(*fi) != NULL)
      if (static_cast<KalmanFitterInfo*>(*fi)->hasReferenceState()) {
        double prevLen = static_cast<KalmanFitterInfo*>(*fi)->getReferenceState()->getForwardSegmentLength();
        double len = prevSegmentLengths[(*fi)->getRep()];
        if (fabs(prevLen + len) > 1E-10 ) {
          std::cerr << "Track::checkConsistency(): segment lengths of reference states don't match" << std::endl;
          std::cerr << prevLen << " + " << len << " = " << prevLen + len << std::endl;
          return false;
        }
        prevSegmentLengths[(*fi)->getRep()] = static_cast<KalmanFitterInfo*>(*fi)->getReferenceState()->getBackwardSegmentLength();
      }
    }

  }

  return true;
}


void Track::trackHasChanged() {

  #ifdef DEBUG
  std::cout << "Track::trackHasChanged \n";
  #endif

  if (fitStatuses_.empty())
    return;

  for (std::map< const AbsTrackRep*, FitStatus* >::const_iterator it=fitStatuses_.begin(); it!=fitStatuses_.end(); ++it) {
    it->second->setHasTrackChanged();
  }
}


void Track::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::Track.
  const bool streamTrackPoints = true; // debugging option
   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::Track thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      {
        std::vector<AbsTrackRep*> &R__stl =  trackReps_;
        R__stl.clear();
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
        if (R__tcl1==0) {
          Error("trackReps_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        R__stl.reserve(R__n);
        for (R__i = 0; R__i < R__n; R__i++) {
          genfit::AbsTrackRep* R__t;
          R__b >> R__t;
          R__stl.push_back(R__t);
        }
      }
      R__b >> cardinalRep_;
      if (streamTrackPoints)
      {
        std::vector<TrackPoint*> &R__stl =  trackPoints_;
        R__stl.clear();
        TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::TrackPoint));
        if (R__tcl1==0) {
          Error("trackPoints_ streamer","Missing the TClass object for genfit::TrackPoint!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        R__stl.reserve(R__n);
        for (R__i = 0; R__i < R__n; R__i++) {
          genfit::TrackPoint* R__t;
          R__t = (genfit::TrackPoint*)R__b.ReadObjectAny(R__tcl1);
          R__t->setTrack(this);
          R__t->fixupRepsForReading();
          R__stl.push_back(R__t);
        }
      }
      {
        std::map<const AbsTrackRep*,FitStatus*> &R__stl =  fitStatuses_;
        R__stl.clear();
        TClass *R__tcl2 = TBuffer::GetClass(typeid(genfit::FitStatus));
        if (R__tcl2==0) {
          Error("fitStatuses_ streamer","Missing the TClass object for genfit::FitStatus!");
          return;
        }
        int R__i, R__n;
        R__b >> R__n;
        for (R__i = 0; R__i < R__n; R__i++) {
          int id;
          R__b >> id;
          genfit::FitStatus* R__t2;
          R__t2 = (genfit::FitStatus*)R__b.ReadObjectAny(R__tcl2);

          R__stl[this->getTrackRep(id)] = R__t2;
        }
      }
      stateSeed_.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      {
        std::vector<AbsTrackRep*> &R__stl =  trackReps_;
        int R__n=(&R__stl) ? int(R__stl.size()) : 0;
        R__b << R__n;
        if(R__n) {
          TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
          if (R__tcl1==0) {
            Error("trackReps_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
            return;
          }
          std::vector<AbsTrackRep*>::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << *R__k;
          }
        }
      }
      R__b << cardinalRep_;
      if (streamTrackPoints)
      {
        std::vector<TrackPoint*> &R__stl =  trackPoints_;
        int R__n=(&R__stl) ? int(R__stl.size()) : 0;
        R__b << R__n;
        if(R__n) {
          std::vector<TrackPoint*>::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
          }
        }
      }
      {
        std::map<const AbsTrackRep*,FitStatus*> &R__stl =  fitStatuses_;
        int R__n=(&R__stl) ? int(R__stl.size()) : 0;
        R__b << R__n;
        if(R__n) {
          TClass *R__tcl1 = TBuffer::GetClass(typeid(genfit::AbsTrackRep));
          if (R__tcl1==0) {
            Error("fitStatuses_ streamer","Missing the TClass object for genfit::AbsTrackRep!");
            return;
          }
          std::map<const AbsTrackRep*,FitStatus*>::iterator R__k;
          for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            int id = this->getIdForRep((*R__k).first);
            R__b << id;
            R__b << ((*R__k).second);
          }
        }
      }
      stateSeed_.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}



} /* End of namespace genfit */
