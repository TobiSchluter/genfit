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

/* This implements the simple Kalman fitter with no reference track
   that uses the stateSeed only until it forgets about it after the
   first few hits.  */

#include "Track.h"
#include "TrackPoint.h"
#include "Exception.h"

#include "KalmanFitter.h"

using namespace genfit;


void KalmanFitter::fitTrack(Track* tr, AbsTrackRep* rep)
{
  std::cout << tr->getNumPoints() << " TrackPoints in this track." << std::endl;
  for (size_t i = 0; i < tr->getNumPoints(); ++i)
    processTrackPoint(tr, tr->getPoint(i), rep, 1);
}


void KalmanFitter::processTrack(Track* tr, AbsTrackRep* rep)
{
  currentState.ResizeTo(tr->getStateSeed());
  currentState = tr->getStateSeed();

  currentState.Print();

  currentCov.ResizeTo(6, 6);
  for (int i = 0; i < 6; i++)
    currentCov(i,i) = 1e4;

  fitTrack(tr, rep);
}

void KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp,
				     AbsTrackRep* rep, int direction)
{
  // Extrapolate to TrackPoint.
  const AbsMeasurement* m = tp->getRawMeasurement(0);
  MeasurementOnPlane mOnPlane = m->constructMeasurementOnPlane(rep);
  const SharedPlanePtr plane = mOnPlane.getPlane();

  MeasuredStateOnPlane state(rep);
  rep->setPosMomCov(&state, currentState, currentCov);
  double extLen = 0;
  try {
    extLen = rep->extrapolateToPlane(&state, plane);
  } catch (genfit::Exception& e) {
    std::cerr << e.what();
    return;
  }

  std::cout << "extrapolated by " << extLen << std::endl;
  const TVectorD& stateVector = state.getState();
  const TVectorD& measurement = mOnPlane.getState();
  stateVector.Print();
  measurement.Print();
  (mOnPlane.getHMatrix() * stateVector).Print();

  // If hit, do Kalman algebra.

  // Store in KalmanInfo asscoiated with the TrackPoint.
}
