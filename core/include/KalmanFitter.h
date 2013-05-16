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

class TrackPoint;

class KalmanFitter : public AbsFitter {
 public:
  KalmanFitter() {}
  ~KalmanFitter() {}

  void fitTrack(Track* tr, AbsTrackRep* rep);

  void processTrack(Track* tr, AbsTrackRep* rep) override;

 private:
  void processTrackPoint(Track* tr, TrackPoint* tp, AbsTrackRep* rep, int direction);
  MeasuredStateOnPlane* currentState;
};

}

#endif //genfit_KalmanFitter_h
