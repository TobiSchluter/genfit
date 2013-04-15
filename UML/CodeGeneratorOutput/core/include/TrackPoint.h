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

/** @addtogroup genfit
 * @{
 */

#ifndef genfit_TrackPoint_h
#define genfit_TrackPoint_h

#include <vector>

#include <TObject.h>

#include "AbsMeasurement.h"
#include "KalmanFitterInfo.h"
#include "MaterialInfo.h"


namespace genfit {

class Track;

class TrackPoint : public TObject {

 public:

  TrackPoint();
  TrackPoint(Track* track);
  TrackPoint(const std::vector< genfit::AbsMeasurement* >& rawMeasurements, Track* track);

  TrackPoint(const TrackPoint&); // copy constructor
  TrackPoint(TrackPoint&&) = default; // move constructor
  TrackPoint& operator=(const TrackPoint&); // assignment operator
  TrackPoint& operator=(TrackPoint&&) = default; // move assignment operator

  ~TrackPoint();


  double getSortingParameter() {return sortingParameter_;}

  Track* const getTrack() const {return track_;}

  std::vector< genfit::AbsMeasurement* >& getRawMeasurements() {return rawMeasurements_;}
  AbsMeasurement* getRawMeasurement(int i);
  unsigned int getNumRawMeasurements() {return rawMeasurements_.size();}
  bool hasRawMeasurements() {return (rawMeasurements_.size() > 0);}

  std::vector< genfit::AbsFitterInfo* >& getFitterInfos() {return fitterInfos_;}
  AbsFitterInfo* getFitterInfo(int i);
  unsigned int getNumFitterInfos() {return fitterInfos_.size();}
  bool hasFitterInfos() {return (fitterInfos_.size() > 0);}

  MaterialInfo* getMaterialInfo() {return material_;}
  bool hasMaterialInfo() {return material_ != nullptr;}


  void setSortingParameter(double sortingParameter) {sortingParameter_ = sortingParameter;}
  void addRawMeasurement(AbsMeasurement* rawMeasurement) {rawMeasurements_.push_back(rawMeasurement);}
  void addFitterInfo(AbsFitterInfo* fitterInfo) {fitterInfos_.push_back(fitterInfo);}
  void deleteFitterInfo(int i);
  void deleteFitterInfo(const AbsTrackRep* rep);
  void setMaterial(MaterialInfo* material);


 private:
  double sortingParameter_;

  /**
   * Pointer to Track where TrackPoint belongs to
   */
  Track* track_; // No ownership

  /** 
   *  Can be more than one, e.g. multiple measurements in the same Si detector, left and right measurements of a wire detector etc.
   * @element-type AbsMeasurement
   */
  std::vector< AbsMeasurement* > rawMeasurements_; // Ownership

  std::vector< AbsFitterInfo* > fitterInfos_; // Ownership

  MaterialInfo* material_; // Ownership


  //ClassDef(TrackPoint,1)

};

} /* End of namespace genfit */

#endif // genfit_TrackPoint_h
