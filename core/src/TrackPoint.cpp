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
  sortingParameter_(0), track_(NULL)//, material_(nullptr)
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

  for (std::vector<AbsMeasurement*>::const_iterator m = rawMeasurements.begin(); m != rawMeasurements.end(); ++m) {
    assert(*m != NULL);
    (*m)->setTrackPoint(this);
    rawMeasurements_.push_back(*m);
  }
}


TrackPoint::TrackPoint(const TrackPoint& rhs) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_)
{
  // clone rawMeasurements
  for (boost::ptr_vector<AbsMeasurement>::const_iterator it = rhs.rawMeasurements_.begin(); it != rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = it->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    for (boost::ptr_vector<AbsFitterInfo>::const_iterator it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      // FIXME: can be solved by clone-ability
      AbsFitterInfo* fi = it2->clone();
      fi->setTrackPoint(this);
      fitterInfos_[(*it).first].push_back(fi);
    }
  }
}

TrackPoint::TrackPoint(const TrackPoint& rhs, const std::map<const AbsTrackRep*, AbsTrackRep*>& map) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_)
{
  // clone rawMeasurements
  for (boost::ptr_vector<AbsMeasurement>::const_iterator it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = it->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    for (boost::ptr_vector<AbsFitterInfo>::const_iterator it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      // FIXME: Again, cloneability can take care of this.
      AbsFitterInfo* fi = it2->clone();
      fi->setRep(map.at(it->first));
      fi->setTrackPoint(this);
      fitterInfos_[map.at(it->first)].push_back(fi);
    }
  }
}


TrackPoint& TrackPoint::operator=(const TrackPoint& rhs) {

  sortingParameter_ = rhs.sortingParameter_;
  track_ = rhs.track_;

  rawMeasurements_.clear();
  fitterInfos_.clear();

  // clone rawMeasurements
  for (boost::ptr_vector<AbsMeasurement>::const_iterator it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = it->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    for (boost::ptr_vector<AbsFitterInfo>::const_iterator it2 = it->second.begin(); it2 != it->second.end();  ++it2 ) {
      // FIXME: cloneability should take care
      AbsFitterInfo* fi = it2->clone();
      fi->setTrackPoint(this);
      fitterInfos_[it->first].push_back(fi);
    }
  }

  return *this;
}


TrackPoint::~TrackPoint() {
  ; // ptr_vector takes care of everything
}


/*void TrackPoint::setMaterial(MaterialInfo* material) {
  if (material_ != nullptr)
    delete material_;
  material_ = material;
}*/


std::vector< genfit::AbsMeasurement* > TrackPoint::getRawMeasurements() {
  // FIXME: Change type?  Why is this needed?
  std::vector<genfit::AbsMeasurement*> vec;
  vec.reserve(rawMeasurements_.size());
  for (size_t i = 0; i < rawMeasurements_.size(); ++i)
    vec.push_back(&rawMeasurements_[i]);
  return vec;

  /*std::vector< genfit::AbsMeasurement* > retVal;
  retVal.reserve(rawMeasurements_.size());

  for (std::vector<AbsMeasurement*>::const_iterator it = rawMeasurements_.begin(); it!=rawMeasurements_.end(); ++it) {
    retVal.push_back(*it);
  }

  return retVal;*/
}


AbsMeasurement* TrackPoint::getRawMeasurement(int i) {
  if (i < 0)
    i += rawMeasurements_.size();

  return &(rawMeasurements_.at(i));
}


std::vector< AbsFitterInfo* > TrackPoint::getFitterInfos() const {
  // FIXME: Should the return type be a ptr_vector?
  std::vector< AbsFitterInfo* > retVal;

  for (std::map<const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      // FIXME: why const? Because of const_iterator.
      retVal.push_back(const_cast<AbsFitterInfo*>(&it->second[i]));
    }
  }

  return retVal;
}

std::vector< AbsFitterInfo* > TrackPoint::getFitterInfos(const AbsTrackRep* rep) const {
  // FIXME: Probably makes more sense to return a reference to the ptr_vector
  std::vector< AbsFitterInfo* > retVal;

  std::map< const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = fitterInfos_.find(rep);

  if (it != fitterInfos_.end()) {
    for (size_t i = 0; i < it->second.size(); ++i) {
      retVal.push_back(const_cast<AbsFitterInfo*>(&it->second[i]));
    }
  }

  return retVal;
}


AbsFitterInfo* TrackPoint::getFitterInfo(const AbsTrackRep* rep, int i) const {
  std::map< const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = fitterInfos_.find(rep);

  if (it != fitterInfos_.end()) {
    if (i < 0)
      i += it->second.size();
    return const_cast<AbsFitterInfo*>(&it->second.at(i));
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
  std::map< const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = fitterInfos_.find(rep);
  if (it == fitterInfos_.end())
    return false;

  return (it->second.size() > 0);
}


void TrackPoint::deleteFitterInfo(AbsTrackRep* rep, int i) {
  std::map< const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::iterator it = fitterInfos_.find(rep);

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
    rawMeasurements_[i].Print();
    std::cout << "............\n";
  }

  for (std::map< const AbsTrackRep*, boost::ptr_vector<AbsFitterInfo> >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {

    for (unsigned int i=0; i<it->second.size(); ++i) {
      std::cout << "FitterInfo Nr. " << i << " for TrackRep " << it->first << "\n";
      it->second[i].Print();
      std::cout << "............\n";
    }

  }


}

} /* End of namespace genfit */
