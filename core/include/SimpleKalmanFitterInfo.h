/* Copyright 2013, Ludwig-Maximilians-Universität München, Technische Universität München
   Authors: Tobias Schlüter & Johannes Rauch

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

#ifndef genfit_SimpleKalmanFitterInfo_h
#define genfit_SimpleKalmanFitterInfo_h

#include <memory>
#include <vector>

#include "TVectorD.h"
#include "TMatrixDSym.h"

#include "AbsFitterInfo.h"

namespace genfit {

class AbsTrackRep;
class TrackPoint;

class SimpleKalmanFitterInfo : public AbsFitterInfo {
  friend class KalmanFitter;   // FIXME: Should be renamed to SimpleKalmanFitter ...

 public:

  SimpleKalmanFitterInfo() {}
  SimpleKalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);

  ~SimpleKalmanFitterInfo() override {};

  //! Deep copy ctor for polymorphic class.
  SimpleKalmanFitterInfo* clone() const override;

  void deleteForwardInfo() override;
  void deleteBackwardInfo() override;
  void deleteReferenceInfo() override { /* exception?  */ };
  void deleteMeasurementInfo() override { /* what should be stored? */};

  MeasurementOnPlane getResidual(bool biased = false, unsigned int iMeasurement = 0) const override;

  virtual void Print(const Option_t* = "") const override;
  
 private:
  std::unique_ptr<MeasuredStateOnPlane> fwPrediction_;
  std::unique_ptr<MeasuredStateOnPlane> bwPrediction_;
  std::vector<MeasurementOnPlane> measurements_;
};

} // namespace genfit  

#endif

// -*- c++ -*-
