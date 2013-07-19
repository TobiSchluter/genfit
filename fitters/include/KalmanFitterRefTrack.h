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

#ifndef genfit_KalmanFitterRefTrack_h
#define genfit_KalmanFitterRefTrack_h

#include "AbsKalmanFitter.h"

namespace genfit {

class KalmanFitterInfo;

class KalmanFitterRefTrack : public AbsKalmanFitter {
 public:
  KalmanFitterRefTrack(unsigned int maxIterations = 4, double deltaPval = 1e-2, double blowUpFactor = 1e3)
    : AbsKalmanFitter(maxIterations, deltaPval, blowUpFactor) {}
  ~KalmanFitterRefTrack() {}

  /**
   * Needs a prepared track!
   */
  void fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction);
  void processTrack(Track* tr, const AbsTrackRep* rep, bool resortHits);
  /**
   * Prepare the track: calc all reference states.
   * If #setSortingParams is true, the extrapolation lengths will be set as sorting parameters
   * of the TrackPoints.
   */
  void prepareTrack(Track* tr, const AbsTrackRep* rep, bool setSortingParams = false);

 private:
  void processTrackPoint(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi, double& chi2, double& ndf, int direction);

  ClassDef(KalmanFitterRefTrack, 1)
};

}  /* End of namespace genfit */
/** @} */

#endif //genfit_KalmanFitterRefTrack_h
