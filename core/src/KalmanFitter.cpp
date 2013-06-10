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
#include "SimpleKalmanFitterInfo.h"

using namespace genfit;


void KalmanFitter::fitTrack(Track* tr, AbsTrackRep* rep, double chi2, size_t ndf, int direction)
{
  chi2 = 0;
  ndf = 0;
  std::cout << tr->getNumPoints() << " TrackPoints in this track." << std::endl;
  for (size_t i = 0; i < tr->getNumPoints(); ++i)
    {
      TrackPoint *tp = 0;
      assert(direction == +1 || direction == -1);
      if (direction == +1)
	tp = tr->getPoint(i);
      else if (direction == -1)
	tp = tr->getPoint(-i-1);
      SimpleKalmanFitterInfo* fi = new SimpleKalmanFitterInfo(tp, rep);
      tp->addFitterInfo(fi);
      processTrackPoint(tr, tp, fi, rep, chi2, ndf, direction);
    }
}


void KalmanFitter::processTrack(Track* tr, AbsTrackRep* rep)
{
  currentState = new MeasuredStateOnPlane(rep);
  TVectorD seed(tr->getStateSeed());
  //seed[0] += 1e-2;  // just so we don't run through the perfectly correct coordinates
  TMatrixDSym cov(6);
  rep->setPosMomCov(currentState, seed, cov);

  currentState->getCov().UnitMatrix();
  currentState->getCov() *= blowUpFactor_;


  double oldChi2FW = 1e6;
  double oldChi2BW = 1e6;
  size_t nIt = 0;
  for(;;) {
    std::cout << "\033[1;21mstate pre" << std::endl;
    currentState->Print();
    std::cout << "\033[0mfitting" << std::endl;
    double chi2FW = 0;
    size_t ndfFW = 0;
    fitTrack(tr, rep, chi2FW, ndfFW, +1);
    std::cout << "\033[1;21mstate post forward" << std::endl;
    currentState->Print();
    std::cout << "\033[0m";

    // Backwards iteration:
    currentState->getCov() *= blowUpFactor_;  // blow up cov
    double chi2BW = 0;
    size_t ndfBW = 0;
    fitTrack(tr, rep, chi2BW, ndfBW, -1);
    std::cout << "\033[1;21mstate post backward" << std::endl;
    currentState->Print();
    std::cout << "\033[0m";

    ++nIt;
    if (nIt > maxIterations_) {
      // FIXME throw exception
      Exception exc("Track fit didn't converge in max iterations.",__LINE__,__FILE__);
      throw exc;
    }
    std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
	      << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
    if (fabs(oldChi2BW - chi2BW) < deltaChi2_)  {
      // Finished
      break;
    }
    else {
      oldChi2BW = chi2BW;
      oldChi2FW = chi2FW;
      currentState->getCov() *= blowUpFactor_;  // blow up cov
    }
  }
  delete currentState;
}

void
KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp, SimpleKalmanFitterInfo* fi,
				AbsTrackRep* rep, double& chi2, size_t& ndf, int direction)
{
  if (!tp->hasRawMeasurements())
    return;

  assert(tp->getNumRawMeasurements() == 1);  // FIXME: should of course support any number

  // Extrapolate to TrackPoint.
  MeasuredStateOnPlane* state = new MeasuredStateOnPlane(*currentState);
  //state.Print();
  if (fi->measurements_.size() == 0) {
    const AbsMeasurement* m = tp->getRawMeasurement(0);
    SharedPlanePtr plane = m->constructPlane(currentState);
    fi->measurements_.push_back(m->constructMeasurementOnPlane(rep, plane));
  }
  const MeasurementOnPlane& mOnPlane = fi->measurements_[0];
  const SharedPlanePtr plane = mOnPlane.getPlane();

  std::cout << "its plane is at R = " << plane->getO().Perp()
	    << " with normal pointing along (" << plane->getNormal().X() << ", " << plane->getNormal().Y() << ", " << plane->getNormal().Z() << ")" << std::endl;

  //state.Print();
  double extLen = 0;
  try {
    extLen = rep->extrapolateToPlane(state, plane);
  } catch (genfit::Exception& e) {
    std::cerr << e.what();
    return;
  }
  std::cout << "extrapolated by " << extLen << std::endl;
  //std::cout << "after extrap: " << std::endl;
  //state.Print();

  // unique_ptr takes care of disposing of the old prediction, takes ownership of state.
  assert(direction == -1 || direction == +1);
  if (direction == +1)
    fi->fwPrediction_ = std::unique_ptr<MeasuredStateOnPlane>(state);
  else if (direction == -1)
    fi->bwPrediction_ = std::unique_ptr<MeasuredStateOnPlane>(state);

  TVectorD stateVector(state->getState());
  TMatrixDSym cov(state->getCov());
  const TVectorD& measurement(mOnPlane.getState());
  const TMatrixDSym& V(mOnPlane.getCov());
  const TMatrixD& H(mOnPlane.getHMatrix());
  std::cout << "State prediction: "; stateVector.Print();
  std::cout << "Cov prediction: "; state->getCov().Print();
  //cov.Print();
  //measurement.Print();

  TVectorD res(measurement - H*stateVector);
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
  std::cout << "Update: "; update.Print();
  //cov.Print();

  stateVector += update;
  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  cov -= covSumInv;
  std::cout << "updated state: "; stateVector.Print();
  std::cout << "updated cov: "; cov.Print();

  TVectorD resNew(measurement - H*stateVector);
  std::cout << "Residual New = (" << resNew(0);
  if (resNew.GetNrows() > 1)
    std::cout << ", " << resNew(1);
  std::cout << ")" << std::endl;

  currentState->setStateCovPlane(stateVector, cov, plane);
  currentState->setAuxInfo(state->getAuxInfo());

  TDecompChol dec(cov);
  TMatrixDSym mist(cov);
  bool status = dec.Invert(mist);
  if (!status) {
    std::cout << "new cov not pos. def." << std::endl;
  }

  // Calculate chi²
  TMatrixDSym HCHt(cov);
  HCHt.Similarity(H);
  HCHt -= V;
  HCHt *= -1;

  TDecompChol decompNew(HCHt);
  TMatrixDSym HCHtInv(decompNew.Invert());

  chi2 += HCHtInv.Similarity(resNew);
  ndf += measurement.GetNrows();
  std::cout << "chi² = " << HCHtInv.Similarity(resNew) << std::endl;
}
