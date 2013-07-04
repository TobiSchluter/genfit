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

#ifndef genfit_AbsMeasurement_h
#define genfit_AbsMeasurement_h

#include "MeasurementOnPlane.h"

#include <TObject.h>


namespace genfit {

class AbsTrackRep;
class TrackPoint;

  /** 
   *  Contains the measurement and covariance in detector coordinates. Detector and hit ids can be used to point back to the original detector hits (clusters etc.).
   */
class AbsMeasurement : public TObject {


 public:

  AbsMeasurement() {};
  AbsMeasurement(int nDims) : rawHitCoords_(nDims), rawHitCov_(nDims) {}
  AbsMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~AbsMeasurement();

  //! Deep copy ctor for polymorphic class.
  virtual AbsMeasurement* clone() const = 0;

  TrackPoint* getTrackPoint() const {return trackPoint_;}
  void setTrackPoint(TrackPoint* tp) {trackPoint_ = tp;}

  const TVectorD& getRawHitCoords() const {return rawHitCoords_;}
  const TMatrixDSym& getRawHitCov() const {return rawHitCov_;}
  int getDetId() const {return detId_;}
  int getHitId() const {return hitId_;}


  /**
   * Construct (virtual) detector plane (use state's TrackRep). It's possible to make corrections to the plane here.
   */
  virtual SharedPlanePtr constructPlane(const StateOnPlane* state) const = 0;

  /**
   * Construct measurement on given plane wrt. given TrackRep. The raw measurement has to be projected onto the plane. It's possible to make corrections to the coordinates here.
   * Usually the vector will contain only one element. But in the case of e.g. WireMeasurements, it will be 2 (left and right).
   */
  virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const AbsTrackRep*, const SharedPlanePtr) const = 0;

  virtual const TMatrixD& getHMatrix(const AbsTrackRep*) const = 0;

  virtual void Print(const Option_t* = "") const;

 protected:

  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
 AbsMeasurement(const AbsMeasurement&);
#ifndef __CINT__
  AbsMeasurement& operator=(const AbsMeasurement&); // default cannot work because TVector and TMatrix = operators don't do resizing
#endif

  TVectorD rawHitCoords_;
  TMatrixDSym rawHitCov_;
  int detId_;
  int hitId_;

  /** 
   *  Pointer to #TrackPoint where the measurement belongs to
   */
  TrackPoint* trackPoint_; // No ownership


  //ClassDef(AbsMeasurement,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsMeasurement_h
