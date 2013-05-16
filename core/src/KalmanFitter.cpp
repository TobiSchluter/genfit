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

#include "TrackPoint.h"

#include "KalmanFitter.h"

using namespace genfit;


void KalmanFitter::fitTrack(Track* tr, AbsTrackRep* rep)
{
  /* How to extrapolate:
  TVector3 pos(0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1),0+gRandom->Gaus(0,0.1));
  TVector3 mom(0,0.5,0.);
  mom *= randomSign();


  genfit::StateOnPlane state(rep);
  rep->setPosMom(&state, pos, mom);

  genfit::SharedPlanePtr origPlane = state.getPlane();
  genfit::StateOnPlane origState(state);

  const TVector3 linePoint(gRandom->Gaus(0,5), gRandom->Gaus(0,5), gRandom->Gaus(0,5));
  const TVector3 lineDirection(gRandom->Gaus(),gRandom->Gaus(),2+gRandom->Gaus());
  const double radius = gRandom->Uniform(10);

  // forth
  try {
    rep->extrapolateToCylinder(&state, radius, linePoint, lineDirection, false);
  }
  catch (genfit::Exception& e) {
    std::cerr << e.what();
  */
}


void KalmanFitter::processTrack(Track* tr, AbsTrackRep* rep)
{
  fitTrack(tr, rep);
}

void KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp,
				     AbsTrackRep* rep, int direction)
{
  // Extrapolate to TrackPoint.
  const AbsMeasurement* m = tp->getRawMeasurement(0);

  // If hit, do Kalman algebra.

  // Store in KalmanInfo asscoiated with the TrackPoint.
}
