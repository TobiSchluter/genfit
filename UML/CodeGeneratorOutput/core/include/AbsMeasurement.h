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

class TrackPoint;

  /** 
   *  Contains the measurement and covariance in detector coordinates. Detector and hit ids can be used to point back to the original detector hits (clusters etc.).
   */
class AbsMeasurement : public TObject {


 public:

  AbsMeasurement();
  AbsMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~AbsMeasurement();

  TrackPoint* const getTrackPoint() const {return trackPoint_;}



  virtual MeasurementOnPlane constructMeasurementOnPlane() = 0;


 protected:
  TVectorD rawHitCoords_;
  TMatrixDSym rawHitCov_;
  int detId_;
  int hitId_;

  /** 
   *  Pointer to #TrackPoint where the measurement belongs to
   */
  TrackPoint* trackPoint_; // No ownership


  ClassDef(AbsMeasurement,1)

};

} /* End of namespace genfit */

#endif // genfit_AbsMeasurement_h
