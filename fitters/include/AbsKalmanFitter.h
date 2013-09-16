/* Copyright 2013, Ludwig-Maximilians Universität München,
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

#ifndef genfit_AbsKalmanFitter_h
#define genfit_AbsKalmanFitter_h

#include "AbsFitter.h"
#include "MeasurementOnPlane.h"

namespace genfit {

class KalmanFitterInfo;

enum eMultipleMeasurementHandling {
  weightedAverage, // weighted average between measurements; used by DAF
  //weightedClosestToReference,
  unweightedClosestToReference,
  //weightedClosestToPrediction,
  unweightedClosestToPrediction
};


class AbsKalmanFitter : public AbsFitter {
 public:
  AbsKalmanFitter(unsigned int maxIterations = 4, double deltaPval = 1e-3, double blowUpFactor = 1e3)
    : maxIterations_(maxIterations), deltaPval_(deltaPval), blowUpFactor_(blowUpFactor), multipleMeasurementHandling_(unweightedClosestToPrediction) {}
  virtual ~AbsKalmanFitter() {;}

  //virtual void fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction) = 0;

  void getChiSquNdf(const Track* tr, const AbsTrackRep* rep, double& bChi2, double& fChi2, double& bNdf,  double& fNdf) const;
  double getChiSqu(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  double getNdf(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  double getRedChiSqu(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  double getPVal(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  eMultipleMeasurementHandling getMultipleMeasurementHandling() const {return multipleMeasurementHandling_;}

  virtual void setMaxIterations(unsigned int n) {maxIterations_ = n;}
  //! How should multiple measurements be handled?
  void setMultipleMeasurementHandling(eMultipleMeasurementHandling mmh) {multipleMeasurementHandling_ = mmh;}

  bool isTrackPrepared(const Track* tr, const AbsTrackRep* rep) const;
  bool isTrackFitted(const Track* tr, const AbsTrackRep* rep) const;

 protected:

  //! get the measurementOnPlane taking the multipleMeasurementHandling_ into account
  const MeasurementOnPlane getMeasurement(const KalmanFitterInfo* fi, int direction) const;

  //! Maximum number of iterations to attempt.  Forward and backward are counted as one iteration.
  unsigned int maxIterations_;
  /**
   * @brief Convergence criterion
   *
   * if track total P-value changes less than this between consecutive iterations, consider the track converged.
   * chi² from the backwards fit is used.
   */
  double deltaPval_;
  //! Blow up the covariance of the forward (backward) fit by this factor before seeding the backward (forward) fit.
  double blowUpFactor_;

  //! How to handle if there are multiple MeasurementsOnPlane
  eMultipleMeasurementHandling multipleMeasurementHandling_;

  ClassDef(AbsKalmanFitter, 1)
};

}
/** @} */

#endif //genfit_AbsKalmanFitter_h
