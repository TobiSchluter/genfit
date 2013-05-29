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

#include <TDecompChol.h>

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
  currentState = new MeasuredStateOnPlane(rep);
  TMatrixDSym cov(6);
  for (int i = 0; i < 6; i++)
    cov(i,i) = 1e2;
  rep->setPosMomCov(currentState, tr->getStateSeed(), cov);

  std::cout << "state pre" << std::endl;
  currentState->Print();
  std::cout << "fitting" << std::endl;
  fitTrack(tr, rep);
  std::cout << "state post" << std::endl;
  currentState->Print();
  delete currentState;
}

void KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp,
				     AbsTrackRep* rep, int direction)
{
  if (!tp->hasRawMeasurements())
    return;

  // Extrapolate to TrackPoint.
  const AbsMeasurement* m = tp->getRawMeasurement(0);
  MeasuredStateOnPlane state(*currentState);
  //state.Print();
  MeasurementOnPlane mOnPlane = m->constructMeasurementOnPlane(rep, *currentState);
  const SharedPlanePtr plane = mOnPlane.getPlane();

  std::cout << "its plane is at R = " << plane->getO().Perp()
	    << " with normal pointing along (" << plane->getNormal().X() << ", " << plane->getNormal().Y() << ", " << plane->getNormal().Z() << ")" << std::endl;

  //state.Print();
  double extLen = 0;
  try {
    extLen = rep->extrapolateToPlane(&state, plane);
  } catch (genfit::Exception& e) {
    std::cerr << e.what();
    return;
  }
  std::cout << "extrapolated by " << extLen << std::endl;
  //std::cout << "after extrap: " << std::endl;
  //state.Print();

  TVectorD stateVector(state.getState());
  TMatrixDSym cov(state.getCov());
  const TVectorD& measurement(mOnPlane.getState());
  const TMatrixDSym& V(mOnPlane.getCov());
  const TMatrixD& H(mOnPlane.getHMatrix());
  stateVector.Print();
  //cov.Print();
  //measurement.Print();

  TVectorD res(measurement - (H * stateVector));
  std::cout << "Residual = (" << res(0);
  if (res.GetNrows() > 1)
    std::cout << ", " << res(1);
  std::cout << ")" << std::endl;
  // If hit, do Kalman algebra.

  // calculate kalman gain ------------------------------
  // calculate covsum (V + HCH^T)
  TMatrixDSym HcovHt(cov);
  HcovHt.Similarity(H);

  TMatrixDSym covSum(V + HcovHt);
  //std::cerr << std::flush << std::endl;
  //std::cout << std::flush;
  //std::cout << "a sum's components:" << std::endl;
  //V.Print();
  //HcovHt.Print();

  TDecompChol decomp(covSum);
  TMatrixDSym covSumInv(decomp.Invert());
  //std::cout << "a matrix and its inverse:" << std::endl;
  //covSum.Print();
  //covSumInv.Print();

  TMatrixD CHt(cov, TMatrixD::kMultTranspose, H);
  TVectorD update = TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res;

  //std::cout << "STATUS:" << std::endl;
  //stateVector.Print();
  update.Print();
  //cov.Print();

  stateVector += update;
  covSumInv.Similarity(CHt);
  cov -= covSumInv;

  TVectorD resNew(measurement - H*stateVector);
  std::cout << "Residual New = (" << resNew(0);
  if (resNew.GetNrows() > 1)
    std::cout << ", " << resNew(1);
  std::cout << ")" << std::endl;

  currentState->setStateCovPlane(stateVector, cov, plane);
  currentState->setAuxInfo(state.getAuxInfo());

  TDecompChol dec(cov);
  TMatrixDSym mist(cov);
  bool status = dec.Invert(mist);
  if (!status)
    {
      std::cout << "new cov not pos. def." << std::endl;
    }

  // Calculate chi²
  TMatrixDSym HCHt(cov);
  HCHt.Similarity(H);
  HCHt -= V;
  HCHt *= -1;

  TDecompChol decompNew(HCHt);
  TMatrixDSym HCHtInv(decompNew.Invert());

  std::cout << "chi² = " << HCHtInv.Similarity(resNew) << std::endl;

  // Store in KalmanInfo asscoiated with the TrackPoint.
}
