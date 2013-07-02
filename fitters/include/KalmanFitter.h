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

#ifndef genfit_KalmanFitter_h
#define genfit_KalmanFitter_h

#include "AbsFitter.h"
#include "MeasuredStateOnPlane.h"

namespace genfit {

class KalmanFitterInfo;
class TrackPoint;

class KalmanFitter : public AbsFitter {
 public:
  KalmanFitter(size_t maxIterations = 4, double deltaChi2 = 1e-3, double blowUpFactor = 1e3)
    : maxIterations_(maxIterations), deltaChi2_(deltaChi2), blowUpFactor_(blowUpFactor) {}
  ~KalmanFitter() {}

  void fitTrack(Track* tr, AbsTrackRep* rep, double chi2, size_t ndf, int direction);

  void processTrack(Track* tr, AbsTrackRep* rep);

 private:
  void processTrackPoint(Track* tr, TrackPoint* tp, KalmanFitterInfo* fi,
			 AbsTrackRep* rep, double& chi2, size_t& ndf, int direction);
  MeasuredStateOnPlane* currentState;

  // Maximum number of iterations to attempt.  Forward and backward
  // are counted as one iteration.
  size_t maxIterations_;
  // Convergence criterion: if track total chi² changes less than this
  // between consecutive iterations, consider the track converged.
  // chi² from the backwards fit is used.
  double deltaChi2_;
  // Blow up the covariance of the forward (backward) fit by this
  // factor before seeding the backward (forward) fit.
  double blowUpFactor_;
};

}
/** @} */

#endif //genfit_KalmanFitter_h
