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

#include "AbsMeasurement.h"
#include "AbsFitterInfo.h"
#include "MaterialInfo.h"

#include <TObject.h>

#include <map>
#include <vector>
#include <memory>

namespace genfit {

class Track;

class TrackPoint : public TObject {

 public:

  TrackPoint();
  TrackPoint(Track* track);

  /**
   * Contructor taking list of measurements. The #setTrackPoint() of each measurement will be called.
   * TrackPoint takes ownership over rawMeasurements.
   */
  TrackPoint(const std::vector< genfit::AbsMeasurement* >& rawMeasurements, Track* track);

  TrackPoint(const TrackPoint&); // copy constructor
  TrackPoint& operator=(const TrackPoint&); // assignment operator

  TrackPoint(const TrackPoint&, const std::map<const AbsTrackRep*, AbsTrackRep*>&); // custom copy constructor where all TrackRep pointers are exchanged according to the map.

  ~TrackPoint();


  double getSortingParameter() const {return sortingParameter_;}

  Track* getTrack() const {return track_;}
  void setTrack(Track* track) {track_ = track;}

  std::vector< genfit::AbsMeasurement* > getRawMeasurements();
  AbsMeasurement* getRawMeasurement(int i = 0) const;
  unsigned int getNumRawMeasurements() const {return rawMeasurements_.size();}
  bool hasRawMeasurements() const {return (! rawMeasurements_.empty());}
  //! Get list of all fitterInfos
  std::vector< AbsFitterInfo* > getFitterInfos() const;
  AbsFitterInfo* getFitterInfo(const AbsTrackRep* rep) const {return fitterInfos_.at(rep);}
  bool hasFitterInfo(const AbsTrackRep* rep) const;

  //MaterialInfo* getMaterialInfo() {return material_;}
  //bool hasMaterialInfo() {return material_ != nullptr;}


  void setSortingParameter(double sortingParameter) {sortingParameter_ = sortingParameter;}
  //! Takes ownership
  void addRawMeasurement(AbsMeasurement* rawMeasurement) {assert(rawMeasurement!=NULL); rawMeasurements_.push_back(rawMeasurement);}
  //! Takes Ownership
  void setFitterInfo(AbsFitterInfo* fitterInfo) {fitterInfos_[fitterInfo->getRep()] = fitterInfo;} // FIXME memory leak
  void deleteFitterInfo(AbsTrackRep* rep) {fitterInfos_.erase(rep);} // FIXME memory leak

  //void setMaterial(MaterialInfo* material);

  void Print(const Option_t* = "") const;

 private:
  double sortingParameter_;

  //! Pointer to Track where TrackPoint belongs to
  Track* track_; // No ownership

  /** 
   *  Can be more than one, e.g. multiple measurements in the same Si detector, left and right measurements of a wire detector etc.
   * @element-type AbsMeasurement
   */
  std::vector<AbsMeasurement*> rawMeasurements_;
  //std::vector<AbsMeasurement*> rawMeasurements_; // Ownership

  std::map< const AbsTrackRep*, AbsFitterInfo* > fitterInfos_; // Ownership over FitterInfos

  //MaterialInfo* material_; // Ownership


  //ClassDef(TrackPoint,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackPoint_h
