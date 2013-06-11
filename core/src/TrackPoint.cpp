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
  sortingParameter_(0), track_(track), rawMeasurements_(rawMeasurements)//, material_(nullptr)
{
  for (AbsMeasurement* measurement : rawMeasurements_) {
    measurement->setTrackPoint(this);
  }
}


TrackPoint::TrackPoint(const TrackPoint& rhs) {
  track_ = rhs.track_;

  for (AbsMeasurement* rhsMeasurement : rhs.rawMeasurements_) {
    AbsMeasurement* tp = rhsMeasurement->clone();
    tp->setTrackPoint(this);
    rawMeasurements_.push_back(tp);
  }

  for (AbsFitterInfo* rhsFitterInfo : rhs.fitterInfos_) {
    AbsFitterInfo* fi = rhsFitterInfo->clone();
    fi->setTrackPoint(this);
    fitterInfos_.push_back(fi);
  }

  //material_ = new MaterialInfo(*(rhs.material_));
}

TrackPoint& TrackPoint::operator=(const TrackPoint& rhs) {
  track_ = rhs.track_;

  for (AbsMeasurement* measurement : rawMeasurements_) {
    if (measurement != nullptr)
      delete measurement;
  }
  rawMeasurements_.clear();

  for (AbsMeasurement* rhsMeasurement : rhs.rawMeasurements_) {
    rawMeasurements_.push_back(rhsMeasurement->clone());
  }

  for (AbsFitterInfo* fitterInfo : fitterInfos_) {
    if (fitterInfo != nullptr)
      delete fitterInfo;
  }
  fitterInfos_.clear();

  for (AbsFitterInfo* rhsFitterInfo : rhs.fitterInfos_) {
    fitterInfos_.push_back(rhsFitterInfo->clone());
  }

  //if (material_ != nullptr)
  //  delete material_;

  return *this;
}


TrackPoint::~TrackPoint() {
  for (AbsMeasurement* measurement : rawMeasurements_) {
    if (measurement != nullptr)
      delete measurement;
  }

  for (AbsFitterInfo* fitterInfo : fitterInfos_) {
    if (fitterInfo != nullptr)
      delete fitterInfo;
  }

  //if (material_ != nullptr)
  //  delete material_;
}


/*void TrackPoint::setMaterial(MaterialInfo* material) {
  if (material_ != nullptr)
    delete material_;
  material_ = material;
}*/


AbsMeasurement* TrackPoint::getRawMeasurement(int i) {
  if (i < 0)
    i += rawMeasurements_.size();

  return rawMeasurements_.at(i);
}


AbsFitterInfo* TrackPoint::getFitterInfo(int i) {
  if (i < 0)
    i += fitterInfos_.size();

  return fitterInfos_.at(i);
}


AbsFitterInfo* TrackPoint::getFitterInfo(const AbsTrackRep* rep, int i) {
  int n(0);
  for (AbsFitterInfo* fitterInfo : fitterInfos_) {
    if (fitterInfo->getRep() == rep) {
      if (i == n) {
        return fitterInfo;
      }
      ++n;
    }
  }

  if (n == 0 || i > n-1) {
    Exception exc("TrackPoint::getFitterInfo ==> cannot get fitter info",__LINE__,__FILE__);
    throw exc;
  }

  if (i<0)
    return getFitterInfo(rep, i+n);

  Exception exc("TrackPoint::getFitterInfo ==> cannot get fitter info",__LINE__,__FILE__);
  throw exc;
}


void TrackPoint::deleteFitterInfo(int i) {
  if (i < 0)
    i += fitterInfos_.size();

  if (fitterInfos_.at(i) != nullptr) {
    delete fitterInfos_[i];
  }
  fitterInfos_.erase(begin(fitterInfos_) + i);
}

void TrackPoint::deleteFitterInfo(const AbsTrackRep* rep) {
  for (unsigned int i=0; i<fitterInfos_.size(); ++i){
    if (fitterInfos_[i]->getRep() == rep) {
      deleteFitterInfo(i);
      --i;
    }
  }
}


void TrackPoint::Print(const Option_t*) const {
  std::cout << "genfit::TrackPoint, belonging to Track " << track_ << "; sorting parameter = " << sortingParameter_ << "\n";
  std::cout << "contains " << rawMeasurements_.size() << " rawMeasurements and " << fitterInfos_.size() << " fitterInfos.\n";

  for (unsigned int i=0; i<rawMeasurements_.size(); ++i) {
    std::cout << "RawMeasurement Nr. " << i << "\n";
    rawMeasurements_[i]->Print();
    std::cout << "............\n";
  }

  for (unsigned int i=0; i<fitterInfos_.size(); ++i) {
    std::cout << "FitterInfo Nr. " << i << "\n";
    fitterInfos_[i]->Print();
    std::cout << "............\n";
  }

}

} /* End of namespace genfit */
