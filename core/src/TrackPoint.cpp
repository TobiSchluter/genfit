/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "TrackPoint.h"

#include "Exception.h"

#include <iostream>


namespace genfit {

TrackPoint::TrackPoint() :
  sortingParameter_(0), track_(nullptr)//, material_(nullptr)
{
  ;
}

TrackPoint::TrackPoint(Track* track) :
  sortingParameter_(0), track_(track)//, material_(nullptr)
{
  ;
}

TrackPoint::TrackPoint(const std::vector< genfit::AbsMeasurement* >& rawMeasurements, Track* track) :
  sortingParameter_(0), track_(track)//, material_(nullptr)
{
  rawMeasurements_.reserve(rawMeasurements.size());

  for (AbsMeasurement* m : rawMeasurements) {
    m->setTrackPoint(this);
    rawMeasurements_.push_back(std::unique_ptr<AbsMeasurement>(m));
  }
}


TrackPoint::TrackPoint(const TrackPoint& rhs) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_)
{
  // clone rawMeasurements
  for (auto it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = it->get()->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (auto it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    for (auto it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      AbsFitterInfo* fi = it2->get()->clone();
      fi->setTrackPoint(this);
      fitterInfos_[it->first].push_back(std::unique_ptr<AbsFitterInfo>(fi));
    }
  }
}

TrackPoint::TrackPoint(const TrackPoint& rhs, const std::map<const AbsTrackRep*, AbsTrackRep*>& map) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_)
{
  // clone rawMeasurements
  for (auto it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = it->get()->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (auto it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    for (auto it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      AbsFitterInfo* fi = it2->get()->clone();
      fi->setRep(map.at(it->first));
      fi->setTrackPoint(this);
      fitterInfos_[map.at(it->first)].push_back(std::unique_ptr<AbsFitterInfo>(fi));
    }
  }
}


TrackPoint& TrackPoint::operator=(const TrackPoint& rhs) {

  sortingParameter_ = rhs.sortingParameter_;
  track_ = rhs.track_;

  rawMeasurements_.clear();
  fitterInfos_.clear();

  // clone rawMeasurements
  for (auto it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = it->get()->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (auto it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    for (auto it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      AbsFitterInfo* fi = it2->get()->clone();
      fi->setTrackPoint(this);
      fitterInfos_[it->first].push_back(std::unique_ptr<AbsFitterInfo>(fi));
    }
  }

  return *this;
}


TrackPoint::~TrackPoint() {
  ; // smart pointers take care of everything
}


/*void TrackPoint::setMaterial(MaterialInfo* material) {
  if (material_ != nullptr)
    delete material_;
  material_ = material;
}*/


std::vector< genfit::AbsMeasurement* > TrackPoint::getRawMeasurements() const {
  std::vector< genfit::AbsMeasurement* > retVal;
  retVal.reserve(rawMeasurements_.size());

  for (auto it = rawMeasurements_.begin(); it!=rawMeasurements_.end(); ++it) {
    retVal.push_back(it->get());
  }

  return retVal;
}


AbsMeasurement* TrackPoint::getRawMeasurement(int i) const {
  if (i < 0)
    i += rawMeasurements_.size();

  return rawMeasurements_.at(i).get();
}


std::vector< AbsFitterInfo* > TrackPoint::getFitterInfos() const {
  std::vector< AbsFitterInfo* > retVal;

  for (auto it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    for (auto it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      retVal.push_back(it2->get());
    }
  }

  return retVal;
}

std::vector< AbsFitterInfo* > TrackPoint::getFitterInfos(const AbsTrackRep* rep) const {
  std::vector< AbsFitterInfo* > retVal;

  auto it = fitterInfos_.find(rep);

  if (it != fitterInfos_.end()) {
    for (auto it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      retVal.push_back(it2->get());
    }
  }

  return retVal;
}


AbsFitterInfo* TrackPoint::getFitterInfo(const AbsTrackRep* rep, int i) const {
  auto it = fitterInfos_.find(rep);

  if (it != fitterInfos_.end()) {
    if (i < 0)
      i += it->second.size();
    return it->second.at(i).get();
  }

  Exception exc("TrackPoint::getFitterInfo ==> no FitterInfo for given TrackRep",__LINE__,__FILE__);
  throw exc;
}


unsigned int TrackPoint::getNumFitterInfos(const AbsTrackRep* rep) const {
  if (!(hasFitterInfos(rep)))
    return 0;

  return fitterInfos_.at(rep).size();
}


bool TrackPoint::hasFitterInfos(const AbsTrackRep* rep) const {
  auto it = fitterInfos_.find(rep);
  if (it == fitterInfos_.end())
    return false;

  return (it->second.size() > 0);
}


void TrackPoint::deleteFitterInfo(AbsTrackRep* rep, int i) {
  auto it = fitterInfos_.find(rep);

  if (it != fitterInfos_.end()) {
    if (i < 0)
      i += it->second.size();

    it->second.erase(it->second.begin()+i);
  }
}

void TrackPoint::deleteFitterInfos(const AbsTrackRep* rep) {
  fitterInfos_.erase(rep);
}


void TrackPoint::Print(const Option_t*) const {
  std::cout << "genfit::TrackPoint, belonging to Track " << track_ << "; sorting parameter = " << sortingParameter_ << "\n";
  std::cout << "contains " << rawMeasurements_.size() << " rawMeasurements and " << getFitterInfos().size() << " fitterInfos for " << fitterInfos_.size() << " TrackReps.\n";

  for (unsigned int i=0; i<rawMeasurements_.size(); ++i) {
    std::cout << "RawMeasurement Nr. " << i << "\n";
    rawMeasurements_[i]->Print();
    std::cout << "............\n";
  }

  for (auto it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {

    for (unsigned int i=0; i<it->second.size(); ++i) {
      std::cout << "FitterInfo Nr. " << i << " for TrackRep " << it->first << "\n";
      it->second[i]->Print();
      std::cout << "............\n";
    }

  }


}

} /* End of namespace genfit */
