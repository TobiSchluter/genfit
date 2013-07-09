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
  for (std::vector<AbsMeasurement*>::const_iterator it = rhs.rawMeasurements_.begin(); it != rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = (*it)->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    AbsFitterInfo* fi = it->second->clone();
    fi->setTrackPoint(this);
    setFitterInfo(fi);
  }
}

TrackPoint::TrackPoint(const TrackPoint& rhs, const std::map<const AbsTrackRep*, AbsTrackRep*>& map) :
  sortingParameter_(rhs.sortingParameter_), track_(rhs.track_)
{
  // clone rawMeasurements
  for (std::vector<AbsMeasurement*>::const_iterator it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = (*it)->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    AbsFitterInfo* fi = it->second->clone();
    fi->setRep(map.at(it->first));
    fi->setTrackPoint(this);
    setFitterInfo(fi);
  }
}


TrackPoint& TrackPoint::operator=(const TrackPoint& rhs) {

  sortingParameter_ = rhs.sortingParameter_;
  track_ = rhs.track_;

  rawMeasurements_.clear();
  fitterInfos_.clear();

  // clone rawMeasurements
  for (std::vector<AbsMeasurement*>::const_iterator it = rhs.rawMeasurements_.begin(); it!=rhs.rawMeasurements_.end(); ++it) {
    AbsMeasurement* tp = (*it)->clone();
    tp->setTrackPoint(this);
    addRawMeasurement(tp);
  }

  // copy fitterInfos
  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = rhs.fitterInfos_.begin(); it != rhs.fitterInfos_.end();  ++it ) {
    AbsFitterInfo* fi = it->second->clone();
    fi->setTrackPoint(this);
    setFitterInfo(fi);
  }

  return *this;
}


TrackPoint::~TrackPoint() {
  // FIXME: We definitely need some smart containers or smart pointers that
  // take care of this, but so far we haven't found a convincing
  // option (2013-07-05).
  
  for (size_t i = 0; i < rawMeasurements_.size(); ++i)
    delete rawMeasurements_[i];

  std::map< const AbsTrackRep*, AbsFitterInfo* >::iterator it;
  for (it = fitterInfos_.begin(); it != fitterInfos_.end(); ++it)
    delete it->second;
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
    vec.push_back(rawMeasurements_[i]);
  return vec;

  /*std::vector< genfit::AbsMeasurement* > retVal;
  retVal.reserve(rawMeasurements_.size());

  for (std::vector<AbsMeasurement*>::const_iterator it = rawMeasurements_.begin(); it!=rawMeasurements_.end(); ++it) {
    retVal.push_back(*it);
  }

  return retVal;*/
}


AbsMeasurement* TrackPoint::getRawMeasurement(int i) const {
  if (i < 0)
    i += rawMeasurements_.size();

  return rawMeasurements_.at(i);
}


std::vector< AbsFitterInfo* > TrackPoint::getFitterInfos() const {
  // FIXME: Should the return type be a ptr_vector?
  std::vector< AbsFitterInfo* > retVal;

  for (std::map<const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    retVal.push_back(it->second);
  }

  return retVal;
}


bool TrackPoint::hasFitterInfo(const AbsTrackRep* rep) const {
  std::map< const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.find(rep);
  if (it == fitterInfos_.end())
    return false;

  return (it->second != NULL);
}


void TrackPoint::Print(const Option_t*) const {
  std::cout << "genfit::TrackPoint, belonging to Track " << track_ << "; sorting parameter = " << sortingParameter_ << "\n";
  std::cout << "contains " << rawMeasurements_.size() << " rawMeasurements and " << getFitterInfos().size() << " fitterInfos for " << fitterInfos_.size() << " TrackReps.\n";

  for (unsigned int i=0; i<rawMeasurements_.size(); ++i) {
    std::cout << "RawMeasurement Nr. " << i << "\n";
    rawMeasurements_[i]->Print();
    std::cout << "............\n";
  }

  for (std::map< const AbsTrackRep*, AbsFitterInfo* >::const_iterator it = fitterInfos_.begin(); it != fitterInfos_.end();  ++it ) {
    std::cout << "FitterInfo for TrackRep " << it->first << "\n";
    it->second->Print();
    std::cout << "............\n";
  }


}

} /* End of namespace genfit */
