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

#ifndef genfit_AbsFitterInfo_h
#define genfit_AbsFitterInfo_h

#include <TObject.h>
#include <TVectorD.h>

#include "MeasurementOnPlane.h"



namespace genfit {

class AbsTrackRep;
class TrackPoint;

  /** 
   *  This class collects all information needed and produced by a specific fitter and is specific to one #GFAbsTrackRep of the #GFTrack.
   */
class AbsFitterInfo : public TObject {

 public:

  AbsFitterInfo();
  AbsFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);

  virtual ~AbsFitterInfo() {};

  //! Deep copy ctor for polymorphic class.
  virtual AbsFitterInfo* clone() const = 0;

  const TrackPoint* getTrackPoint() const {return trackPoint_;}
  const AbsTrackRep* getRep() const {return rep_;}

  void setTrackPoint(const TrackPoint *tp) {trackPoint_ = tp;}
  void setStatusFlag(int statusFlag) {statusFlag_ = statusFlag;}
  virtual void setRep(const AbsTrackRep* rep) {rep_ = rep;}

  virtual void deleteForwardInfo() = 0;
  virtual void deleteBackwardInfo() = 0;
  virtual void deleteReferenceInfo() = 0;
  virtual void deleteMeasurementInfo() = 0;

  virtual MeasuredStateOnPlane getFittedState(bool biased = false) const = 0;
  virtual MeasurementOnPlane getResidual(bool biased = false, unsigned int iMeasurement = 0) const = 0;
  int getStatusFlag() const {return statusFlag_;}

  virtual void Print(const Option_t* = "") const {;}

  virtual bool checkConsistency() const = 0;

 protected:

  /** Pointer to #TrackPoint where the FitterInfo belongs to
   */
  const TrackPoint* trackPoint_; // No ownership

  /** Pointer to #TrackRep with respect to which the FitterInfo is defined
   */
  const AbsTrackRep* rep_; // No ownership

  int statusFlag_; // 0 = ok

  //ClassDef(AbsFitterInfo,1)

 private:
  AbsFitterInfo(const AbsFitterInfo&); // copy constructor
  AbsFitterInfo& operator=(const AbsFitterInfo&); // assignment operator
};

// Needed for boost cloneability:
inline AbsFitterInfo* new_clone( const AbsFitterInfo & a)
{
  return a.clone();
}


} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsFitterInfo_h
