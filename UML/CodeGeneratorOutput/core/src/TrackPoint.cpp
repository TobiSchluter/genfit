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

namespace genfit {

TrackPoint::TrackPoint() :
  sortingParameter_(0), track_(nullptr), material_(nullptr)
{
  ;
}

TrackPoint::TrackPoint(Track* track) :
  sortingParameter_(0), track_(track), material_(nullptr)
{
  ;
}

TrackPoint::TrackPoint(const std::vector< genfit::AbsMeasurement* >& rawMeasurements, Track* track) :
  sortingParameter_(0), track_(track), rawMeasurements_(rawMeasurements), material_(nullptr)
{
  ;
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

  if (material_ != nullptr)
    delete material_;

}


void TrackPoint::setMaterial(MaterialInfo* material) {
  if (material_ != nullptr)
    delete material_;
  material_ = material;
}

} /* End of namespace genfit */
