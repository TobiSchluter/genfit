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

#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"
#include "Exception.h"

#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"

#define DEBUG

using namespace genfit;


void KalmanFitterRefTrack::fitTrack(Track* tr, AbsTrackRep* rep, double chi2, size_t ndf, int direction)
{
  chi2 = 0;
  ndf = 0;
  KalmanFitterInfo* prevFi(nullptr);
  std::cout << tr->getNumPoints() << " TrackPoints with measurements in this track." << std::endl;
  for (size_t i = 0; i < tr->getNumPointsWithMeasurement(); ++i) {
      TrackPoint *tp = 0;
      assert(direction == +1 || direction == -1);
      if (direction == +1)
        tp = tr->getPointWithMeasurement(i);
      else if (direction == -1)
        tp = tr->getPointWithMeasurement(-i-1);

      KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(-1));
      processTrackPoint(fi, prevFi, chi2, ndf, direction);
      prevFi = fi;
  }
}


void KalmanFitterRefTrack::processTrack(Track* tr, AbsTrackRep* rep)
{

  // TODO: try catch block, what if fit fails?

  for (unsigned int i=0; i<maxIterations_; ++i) {

#ifdef DEBUG
    std::cout << " KalmanFitterRefTrack::processTrack, iteration nr. " << i << "\n";
#endif

    // prepare
    prepareTrack(tr, rep);

#ifdef DEBUG
    std::cout << "Prepared Track:"; tr->Print();
#endif

    // fit forward
    double oldChi2FW = 1e6;
    double oldChi2BW = 1e6;
    double chi2FW = 0;
    size_t ndfFW = 0;

    fitTrack(tr, rep, chi2FW, ndfFW, +1);

    // fit backward
    KalmanFitterInfo* lastInfo = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(-1)->getFitterInfo(-1));
    lastInfo->setBackwardPrediction(new MeasuredStateOnPlane(*(lastInfo->getForwardUpdate())));
    lastInfo->getBackwardPrediction()->getCov() *= blowUpFactor_;  // blow up cov
    double chi2BW = 0;
    size_t ndfBW = 0;

    fitTrack(tr, rep, chi2BW, ndfBW, -1);


#ifdef DEBUG
    std::cout << "Track after fit:"; tr->Print();
#endif


    std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
        << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
    if (fabs(oldChi2BW - chi2BW) < deltaChi2_)  {
      // Finished
      break;
    }
    else {
      oldChi2BW = chi2BW;
      oldChi2FW = chi2FW;
    }

  }

}


void KalmanFitterRefTrack::prepareTrack(Track* tr, AbsTrackRep* rep) {

  std::unique_ptr<MeasuredStateOnPlane> seedState;

  // get seed state from previous fit if there is one
  if (tr->getPointWithMeasurement(0)->hasFitterInfos()) {
    KalmanFitterInfo* fitterInfo(nullptr);

    // get the last fitter info with the correct TrackRep and see if it has the right type
    std::vector< genfit::AbsFitterInfo* >& fitterInfos = tr->getPointWithMeasurement(0)->getFitterInfos();
    for (int i=fitterInfos.size()-1; i>-1; --i) {
      if (fitterInfos[i]->getRep() == rep)
        fitterInfo = dynamic_cast<KalmanFitterInfo*>(fitterInfos[i]);
      if (fitterInfo) {
        if (fitterInfo->hasBackwardUpdate()) {
          seedState = std::unique_ptr<MeasuredStateOnPlane>(new MeasuredStateOnPlane(*(fitterInfo->getBackwardUpdate())));
          break;
        }
      }
    }
  }
  // else create seed state from seed info of track
  if (!seedState) {
    seedState = std::unique_ptr<MeasuredStateOnPlane>(new MeasuredStateOnPlane(rep));
    TMatrixDSym cov(rep->getDim());
    cov.UnitMatrix();
    cov *= blowUpFactor_;
    rep->setPosMom(&*seedState, tr->getStateSeed());
    seedState->setCov(cov);
  }

#ifdef DEBUG
  std::cout << "seed state  "; seedState->Print();
#endif

  TMatrixD FTransportMatrix(rep->getDim(), rep->getDim());
  FTransportMatrix.UnitMatrix();
  TMatrixD BTransportMatrix(rep->getDim(), rep->getDim());

  TMatrixDSym FNoiseMatrix(rep->getDim());
  TMatrixDSym BNoiseMatrix(rep->getDim());

  KalmanFitterInfo* prevFitterInfo(nullptr);

  for (unsigned int i=0; i<tr->getNumPoints(); ++i){
    TrackPoint* trackPoint = tr->getPoint(i);

    // check if we have a measurement
    if (!trackPoint->hasRawMeasurements())
      continue;

    // create new fitterInfo and ReferenceState
    KalmanFitterInfo* fitterInfo = new KalmanFitterInfo(trackPoint, rep);
    trackPoint->addFitterInfo(fitterInfo);
    ReferenceStateOnPlane* refState = new ReferenceStateOnPlane(*seedState);
    fitterInfo->setReferenceState(refState);

    // Construct plane
    SharedPlanePtr plane = trackPoint->getRawMeasurement(0)->constructPlane(&*seedState);

    // do extrapolation and set reference state infos
    double segmentLen = rep->extrapolateToPlane(&*seedState, plane);
    if (i>0) rep->getForwardJacobianAndNoise(FTransportMatrix, FNoiseMatrix);
    rep->getBackwardJacobianAndNoise(BTransportMatrix, BNoiseMatrix);

    if (i==0) { // if we are at first measurement and seed state is defined somewhere else, still set forward info to default
      segmentLen = 0;
    }

    refState->setForwardSegmentLength(segmentLen);
    refState->setForwardTransportMatrix(FTransportMatrix);
    refState->setForwardNoiseMatrix(FNoiseMatrix);

    if (prevFitterInfo != nullptr) {
      ReferenceStateOnPlane* prevRefState =  prevFitterInfo->getReferenceState();
      prevRefState->setBackwardSegmentLength(-segmentLen);
      prevRefState->setBackwardTransportMatrix(BTransportMatrix);
      prevRefState->setBackwardNoiseMatrix(BNoiseMatrix);
    }

    prevFitterInfo = fitterInfo;


    // set seed as prediction if at first measurement
    if (i==0) {
      fitterInfo->setForwardPrediction(new MeasuredStateOnPlane(*seedState));
    }

    // get MeasurementsInPlane
    for (AbsMeasurement* measurement : trackPoint->getRawMeasurements()) {
      fitterInfo->addMeasurementOnPlane(new MeasurementOnPlane(measurement->constructMeasurementOnPlane(rep, plane)));
    }

  }

  // set backward info for last reference state
  ReferenceStateOnPlane* prevRefState =  prevFitterInfo->getReferenceState();
  prevRefState->setBackwardSegmentLength(0);
  BTransportMatrix.UnitMatrix();
  prevRefState->setBackwardTransportMatrix(BTransportMatrix);
  BNoiseMatrix.Zero();
  prevRefState->setBackwardNoiseMatrix(BNoiseMatrix);
}


void
KalmanFitterRefTrack::processTrackPoint(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi, double& chi2, size_t& ndf, int direction)
{

  const MeasurementOnPlane& m = fi->getAvgWeightedMeasurementOnPlane();

  TVectorD dp(m.getHMatrix().GetNcols()); // \delta p_{k|k-1}
  TMatrixDSym C(m.getHMatrix().GetNcols()); // C_{k|k-1}

  // predict
  if (prevFi != nullptr) {
    const TMatrixD& F = fi->getReferenceState()->getTransportMatrix(direction); // Transport matrix
    const TMatrixDSym& N = fi->getReferenceState()->getNoiseMatrix(direction); // Noise matrix
    dp = F * (prevFi->getUpdate(direction)->getState() - prevFi->getReferenceState()->getState());
    C = prevFi->getUpdate(direction)->getCov();
    C.Similarity(F);
    C += N;
    fi->setPrediction(new MeasuredStateOnPlane(dp + fi->getReferenceState()->getState(), C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep()), direction);
#ifdef DEBUG
    std::cout << "\033[31m";
    std::cout << "p_{k-1,r} "; prevFi->getReferenceState()->getState().Print();
#endif
  }
  else {
    dp = fi->getPrediction(direction)->getState() - fi->getReferenceState()->getState();
    C = fi->getPrediction(direction)->getCov();
#ifdef DEBUG
    std::cout << "\033[31m";
    std::cout << "p_{k,r} "; fi->getReferenceState()->getState().Print();
#endif
  }

#ifdef DEBUG
  std::cout << "Δp_{k|k-1} "; dp.Print();
  std::cout << " p_{k|k-1} "; fi->getPrediction(direction)->getState().Print();
  std::cout << " C_{k|k-1} "; C.Print();
  std::cout << "\033[0m";
#endif

  // update
  const TVectorD& dm = m.getState() - (m.getHMatrix() * fi->getReferenceState()->getState()); // \delta m_k = m_k - H_k p_{k,r}

  TMatrixDSym covSumInv(C); // (V_k + H_k C_{k|k-1} H_k^T)^(-1)
  covSumInv.Similarity(m.getHMatrix());
  covSumInv += m.getCov();
  tools::invertMatrix(covSumInv);

  TMatrixD CHt(C, TMatrixD::kMultTranspose, m.getHMatrix());

  const TVectorD& res = dm - m.getHMatrix()*dp; //
  //TVectorD res(dm); res.Zero(); // FIXME this is only a test!
#ifdef DEBUG
  std::cout << "\033[34m";
  std::cout << "m   "; m.getState().Print();
  std::cout << "res "; res.Print();
  std::cout << "\033[0m";
#endif
  TVectorD updated(TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
  updated += dp;
  updated += fi->getReferenceState()->getState();

  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  C -= covSumInv; // updated Cov

#ifdef DEBUG
  std::cout << "\033[32m";
  std::cout << "Δp_{k|k} "; (updated - fi->getReferenceState()->getState()).Print();
  std::cout << " p_{k|k} "; updated.Print();
  std::cout << " C_{k|k} "; C.Print();
  std::cout << "\033[0m";
#endif

  // Calculate chi²
  TMatrixDSym Rinv(C);
  Rinv.Similarity(m.getHMatrix());
  Rinv -= m.getCov();
  Rinv *= -1;
  tools::invertMatrix(Rinv);

  TVectorD resNew(m.getState() - m.getHMatrix()*updated);

  double chi2inc = Rinv.Similarity(resNew);
  chi2 += chi2inc;

  double ndfInc = m.getState().GetNrows() * m.getWeight();
  ndf += ndfInc;

#ifdef DEBUG
  std::cout << " chi² inc " << chi2inc << "\n";
  std::cout << " ndf inc  " << ndfInc << "\n";
#endif


  fi->setUpdate(new KalmanFittedStateOnPlane(updated, C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), chi2inc, ndfInc), direction);

}
