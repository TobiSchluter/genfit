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

#include "KalmanFitter.h"

#include "Exception.h"
#include "KalmanFitterInfo.h"
#include "Track.h"
#include "TrackPoint.h"
#include "Tools.h"



//#define DEBUG

using namespace genfit;


void KalmanFitter::fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{

  if (multipleMeasurementHandling_ == unweightedClosestToReference) {
    Exception exc("KalmanFitter::fitTrack ==> cannot use unweightedClosestToReference as multiple measurement handling.",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  chi2 = 0;
  ndf = 0;
#ifdef DEBUG
  std::cout << tr->getNumPoints() << " TrackPoints in this track." << std::endl;
#endif
  for (size_t i = 0; i < tr->getNumPoints(); ++i) {
    TrackPoint *tp = 0;
    KalmanFitterInfo* fi;
    assert(direction == +1 || direction == -1);
    if (direction == +1) {
      tp = tr->getPoint(i);
    }
    else {
      tp = tr->getPoint(-i-1);
    }

    if (tp->getNumFitterInfos(rep) == 0) {
      fi = new KalmanFitterInfo(tp, rep);
      tp->addFitterInfo(fi);
    }
    else
      fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));

    try {
      processTrackPoint(tr, tp, fi, rep, chi2, ndf, direction);
    }
    catch (genfit::Exception& e) {
      fi->setStatusFlag(1);
      std::cerr << e.what();
      return;
    }

  }
}


void KalmanFitter::processTrack(Track* tr, const AbsTrackRep* rep)
{
  currentState = new MeasuredStateOnPlane(rep);
  TVectorD seed(tr->getStateSeed());
  //seed[0] += 1e-2;  // just so we don't run through the perfectly correct coordinates
  TMatrixDSym cov(6);
  rep->setPosMomCov(currentState, seed, cov);

  currentState->getCov().UnitMatrix();
  currentState->getCov() *= blowUpFactor_;

#ifdef DEBUG
  double oldChi2FW = 1e6;
#endif
  double oldChi2BW = 1e6;
  size_t nIt = 0;
  for(;;) {
#ifdef DEBUG
    std::cout << "\033[1;21mstate pre" << std::endl;
    currentState->Print();
    std::cout << "\033[0mfitting" << std::endl;
#endif
    double chi2FW = 0;
    double ndfFW = 0;
    fitTrack(tr, rep, chi2FW, ndfFW, +1);
#ifdef DEBUG
    std::cout << "\033[1;21mstate post forward" << std::endl;
    currentState->Print();
    std::cout << "\033[0m";
#endif

    // Backwards iteration:
    currentState->getCov() *= blowUpFactor_;  // blow up cov
    double chi2BW = 0;
    double ndfBW = 0;
    fitTrack(tr, rep, chi2BW, ndfBW, -1);
#ifdef DEBUG
    std::cout << "\033[1;21mstate post backward" << std::endl;
    currentState->Print();
    std::cout << "\033[0m";

    std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
	      << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
#endif

    if (fabs(oldChi2BW - chi2BW) < deltaChi2_)  {
      // Finished
      break;
    }
    else {
      oldChi2BW = chi2BW;
#ifdef DEBUG
      oldChi2FW = chi2FW;
#endif
      currentState->getCov() *= blowUpFactor_;  // blow up cov
    }

    if (++nIt > maxIterations_) {
      break;
      // FIXME
      //delete currentState;
      //Exception exc("Track fit didn't converge in max iterations.",__LINE__,__FILE__);
      //throw exc;
    }
  }
  delete currentState;
}

void
KalmanFitter::processTrackPoint(Track* tr, TrackPoint* tp, KalmanFitterInfo* fi,
    const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{
  if (!tp->hasRawMeasurements())
    return;

  // Extrapolate to TrackPoint.
  MeasuredStateOnPlane* state = new MeasuredStateOnPlane(*currentState);
  //state.Print();

  // construct measurementsOnPlane if it has not yet been done
  if (fi->getNumMeasurements() == 0) {
    std::vector< genfit::AbsMeasurement* > rawMeasurements =  tp->getRawMeasurements();
    // construct plane with first measurement
    SharedPlanePtr plane = rawMeasurements[0]->constructPlane(currentState);
    for (std::vector< genfit::AbsMeasurement* >::iterator it = rawMeasurements.begin(); it != rawMeasurements.end(); ++it) {
      fi->addMeasurementsOnPlane((*it)->constructMeasurementsOnPlane(rep, plane));
    }
  }
  const SharedPlanePtr plane = fi->getPlane();

#ifdef DEBUG
  std::cout << "its plane is at R = " << plane->getO().Perp()
	    << " with normal pointing along (" << plane->getNormal().X() << ", " << plane->getNormal().Y() << ", " << plane->getNormal().Z() << ")" << std::endl;
#endif

  //state.Print();

#ifdef DEBUG
  double extLen = rep->extrapolateToPlane(state, plane);
  std::cout << "extrapolated by " << extLen << std::endl;
#else
  rep->extrapolateToPlane(state, plane);
#endif
  //std::cout << "after extrap: " << std::endl;
  //state.Print();

  // unique_ptr takes care of disposing of the old prediction, takes ownership of state.
  assert(direction == -1 || direction == +1);
  fi->setPrediction(state, direction);

  TVectorD stateVector(state->getState());
  TMatrixDSym cov(state->getCov());
  const MeasurementOnPlane& mOnPlane = getMeasurement(fi, direction);
  const TVectorD& measurement(mOnPlane.getState());
  const TMatrixDSym& V(mOnPlane.getCov());
  const TMatrixD& H(mOnPlane.getHMatrix());
#ifdef DEBUG
  std::cout << "State prediction: "; stateVector.Print();
  std::cout << "Cov prediction: "; state->getCov().Print();
  //cov.Print();
  //measurement.Print();
#endif

  TVectorD res(measurement - H*stateVector);
#ifdef DEBUG
  std::cout << "Residual = (" << res(0);
  if (res.GetNrows() > 1)
    std::cout << ", " << res(1);
  std::cout << ")" << std::endl;
#endif
  // If hit, do Kalman algebra.

  // calculate kalman gain ------------------------------
  // calculate covsum (V + HCH^T)
  TMatrixDSym HcovHt(cov);
  HcovHt.Similarity(H);

  TMatrixDSym covSumInv(V + HcovHt);
  tools::invertMatrix(covSumInv);

  TMatrixD CHt(cov, TMatrixD::kMultTranspose, H);
  TVectorD update = TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res;

#ifdef DEBUG
  //std::cout << "STATUS:" << std::endl;
  //stateVector.Print();
  std::cout << "Update: "; update.Print();
  //cov.Print();
#endif

  stateVector += update;
  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  cov -= covSumInv;
#ifdef DEBUG
  std::cout << "updated state: "; stateVector.Print();
  std::cout << "updated cov: "; cov.Print();
#endif

  TVectorD resNew(measurement - H*stateVector);
#ifdef DEBUG
  std::cout << "Residual New = (" << resNew(0);

  if (resNew.GetNrows() > 1)
    std::cout << ", " << resNew(1);
  std::cout << ")" << std::endl;
#endif

  currentState->setStateCovPlane(stateVector, cov, plane);
  currentState->setAuxInfo(state->getAuxInfo());

  /*TDecompChol dec(cov);
  TMatrixDSym mist(cov);
  bool status = dec.Invert(mist);
  if (!status) {
#ifdef DEBUG
    std::cout << "new cov not pos. def." << std::endl;
#endif
  }*/

  // Calculate chi²
  TMatrixDSym HCHt(cov);
  HCHt.Similarity(H);
  HCHt -= V;
  HCHt *= -1;

  tools::invertMatrix(HCHt);

  double chi2inc = HCHt.Similarity(resNew);
  chi2 += chi2inc;

  double ndfInc = measurement.GetNrows();
  ndf += ndfInc;
#ifdef DEBUG
  std::cout << "chi² increment = " << chi2inc << std::endl;
#endif

  // set update
  KalmanFittedStateOnPlane* updatedSOP = new KalmanFittedStateOnPlane(*currentState, chi2inc, ndfInc);
  fi->setUpdate(updatedSOP, direction);
}
