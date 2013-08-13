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
class TrackPoint;

class KalmanFitterRefTrack : public AbsKalmanFitter {
 public:
  KalmanFitterRefTrack(unsigned int maxIterations = 4, double deltaPval = 1e-3, double blowUpFactor = 1e3)
    : AbsKalmanFitter(maxIterations, deltaPval, blowUpFactor), refitAll_(false), deltaChi2Ref_(1) {}
  ~KalmanFitterRefTrack() {cleanSOPsToDestruct();}

  /**
   * Needs a prepared track! Return last TrackPoint that has been processed.
   */
  TrackPoint* fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction);
  void processTrack(Track* tr, const AbsTrackRep* rep, bool resortHits);
  /**
   * Prepare the track: calc all reference states.
   * If #setSortingParams is true, the extrapolation lengths will be set as sorting parameters
   * of the TrackPoints.
   * Returns if the track has been changed.
   */
  bool prepareTrack(Track* tr, const AbsTrackRep* rep, bool setSortingParams = false);


  void setRefitAll(bool refit = true) {refitAll_ = refit;}

  /**
   * When will the reference track be updated?
   * If (smoothedState - referenceState) * smoothedCov^(-1) * (smoothedState - referenceState)^T >= deltaChi2Ref_.
   */
  void setDeltaChi2Ref(double dChi2) {deltaChi2Ref_ = dChi2;}

 private:
  void processTrackPoint(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi, double& chi2, double& ndf, int direction);

  /**
   * Remove referenceStates if they are too far from smoothed states.
   * Does NOT remove forward and backward info, but returns from/to where they have to be removed later
   * Return if anything has changed.
   */
  bool removeOutdated(Track* tr, const AbsTrackRep* rep, int& notChangedUntil, int& notChangedFrom) const;

  /**
   * If refitAll_, remove all information.
   */
  void removeForwardBackwardInfo(Track* tr, const AbsTrackRep* rep, int notChangedUntil, int notChangedFrom) const;

  void cleanSOPsToDestruct();

  bool refitAll_; // always refit all points or only if reference states have changed
  double deltaChi2Ref_; // reference track update cut

  std::vector<const StateOnPlane*> SOPsToDestruct_; //! ownership, "garbage collector" helper

 public:
  ClassDef(KalmanFitterRefTrack, 1)

};

}  /* End of namespace genfit */
/** @} */

#endif //genfit_KalmanFitterRefTrack_h
