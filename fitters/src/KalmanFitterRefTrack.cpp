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
#include "KalmanFitStatus.h"

#include <Math/ProbFunc.h>


//#define DEBUG


using namespace genfit;


void KalmanFitterRefTrack::fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{

  if (!isTrackPrepared(tr, rep)) {
    Exception exc("KalmanFitterRefTrack::fitTrack ==> track is not properly prepared.",__LINE__,__FILE__);
    throw exc;
  }

  chi2 = 0;
  ndf = 0;
  KalmanFitterInfo* prevFi(NULL);

#ifdef DEBUG
  std::cout << tr->getNumPoints() << " TrackPoints with measurements in this track." << std::endl;
#endif

  for (size_t i = 0; i < tr->getNumPointsWithMeasurement(); ++i) {
      TrackPoint *tp = 0;
      assert(direction == +1 || direction == -1);
      if (direction == +1)
        tp = tr->getPointWithMeasurement(i);
      else if (direction == -1)
        tp = tr->getPointWithMeasurement(-i-1);

      KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));
      processTrackPoint(fi, prevFi, chi2, ndf, direction);

      prevFi = fi;
  }
}


void KalmanFitterRefTrack::processTrack(Track* tr, const AbsTrackRep* rep)
{

#ifdef DEBUG
  double oldChi2FW = 1e6;
  double oldPvalFW = 0.;
  double oldChi2BW = 1e6;
#endif
  double oldPvalBW = 0.;
  double chi2FW(0), ndfFW(0);
  double chi2BW(0), ndfBW(0);

  KalmanFitStatus* status = new KalmanFitStatus();
  tr->setFitStatus(status, rep);

  status->setIsFittedWithReferenceTrack(true);

  unsigned int nIt=0;
  for (;;) {

    try {
      #ifdef DEBUG
      std::cout << " KalmanFitterRefTrack::processTrack, iteration nr. " << i << "\n";
      #endif

      // prepare
      prepareTrack(tr, rep);

      #ifdef DEBUG
      std::cout << "Prepared Track:"; tr->Print();
      #endif

      // fit forward
      #ifdef DEBUG
      std::cout << "forward fit\n";
      #endif
      fitTrack(tr, rep, chi2FW, ndfFW, +1);

      // fit backward
      #ifdef DEBUG
      std::cout << "backward fit\n";
      #endif
      KalmanFitterInfo* lastInfo = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(-1)->getFitterInfo(rep));
      lastInfo->setBackwardPrediction(new MeasuredStateOnPlane(*(lastInfo->getForwardUpdate())));
      lastInfo->getBackwardPrediction()->getCov() *= blowUpFactor_;  // blow up cov

      fitTrack(tr, rep, chi2BW, ndfBW, -1);

      ++nIt;


      #ifdef DEBUG
      std::cout << "Track after fit:"; tr->Print();


      std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
          << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
      #endif

      // See if p-value only changed little.  If the initial
      // parameters are very far off, initial chi^2 and the chi^2
      // after the first iteration will be both very close to zero, so
      // we need to force at least two iterations here.  Convergence
      // doesn't make much sense before running twice anyway.
      double PvalBW = ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW);
      #ifdef DEBUG
      double PvalFW = ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW);
      #endif
      if (nIt > 1 && fabs(oldPvalBW - PvalBW) < deltaPval_)  {
        // Finished
        status->setIsFitConverged();
        break;
      }
      else {
        oldPvalBW = PvalBW;
        #ifdef DEBUG
        oldChi2BW = chi2BW;
        oldChi2FW = chi2FW;
        oldPvalFW = PvalFW;
        #endif
      }

      if (nIt >= maxIterations_) {
        break;
      }
    }
    catch(Exception& e) {
      std::cerr << e.what();
      status->setIsFitted(false);
      status->setIsFitConverged(false);
      return;
    }

  }


  // check
  assert(tr->checkConsistency());

  status->setIsFitted();
  status->setCharge(rep->getCharge(static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
  status->setNumIterations(nIt);
  status->setForwardChiSqu(chi2FW);
  status->setBackwardChiSqu(chi2BW);
  status->setForwardNdf(ndfFW);
  status->setBackwardNdf(ndfBW);

}


void KalmanFitterRefTrack::prepareTrack(Track* tr, const AbsTrackRep* rep) {

  std::auto_ptr<MeasuredStateOnPlane> seedState;

  // get seed state from previous fit if there is one
  if (tr->getPointWithMeasurement(0)->hasFitterInfo(rep)) {

    // get the last fitter info with the correct TrackRep and see if it has the right type
    AbsFitterInfo* fi = tr->getPointWithMeasurement(0)->getFitterInfo(rep);
    KalmanFitterInfo* fitterInfo = dynamic_cast<KalmanFitterInfo*>(fi);
    if (fitterInfo) {
      if (fitterInfo->hasBackwardUpdate()) {
        seedState = std::auto_ptr<MeasuredStateOnPlane>(new MeasuredStateOnPlane(*(fitterInfo->getBackwardUpdate())));
        seedState->getCov() *= blowUpFactor_;
      }
    }
  }
  // else create seed state from seed info of track
  if (seedState.get() == NULL) {
    seedState = std::auto_ptr<MeasuredStateOnPlane>(new MeasuredStateOnPlane(rep));
    rep->setPosMom(&*seedState, tr->getStateSeed());
    TMatrixDSym cov(rep->getDim());
    cov.UnitMatrix();
    //cov *= blowUpFactor_; // FIXME find good start values
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

  KalmanFitterInfo* prevFitterInfo(NULL);

  for (unsigned int i=0; i<tr->getNumPoints(); ++i){
    TrackPoint* trackPoint = tr->getPoint(i);

    // check if we have a measurement
    if (!trackPoint->hasRawMeasurements())
      continue;

    // create new fitterInfo and ReferenceState
    KalmanFitterInfo* fitterInfo;
    AbsFitterInfo* absFitterInfo;
    std::vector<double> oldWeights;
    bool KalmanFiAvailable = false;
    if (trackPoint->hasFitterInfo(rep)) {
      absFitterInfo = trackPoint->getFitterInfo(rep);
      if (dynamic_cast<KalmanFitterInfo*>(absFitterInfo) != NULL) {
        KalmanFiAvailable = true;
      }
    }
    if (KalmanFiAvailable) {
      fitterInfo = static_cast<KalmanFitterInfo*>(absFitterInfo);
      oldWeights = fitterInfo->getWeights();
      fitterInfo->deleteForwardInfo();
      fitterInfo->deleteBackwardInfo();
      fitterInfo->deleteReferenceInfo();
      fitterInfo->deleteMeasurementInfo();
    }
    else {
      fitterInfo = new KalmanFitterInfo(trackPoint, rep);
      trackPoint->setFitterInfo(fitterInfo);
    }


    // Construct plane
    SharedPlanePtr plane = trackPoint->getRawMeasurement(0)->constructPlane(&*seedState);

    // do extrapolation and set reference state infos
    double segmentLen = rep->extrapolateToPlane(&*seedState, plane);
    if (i>0) rep->getForwardJacobianAndNoise(FTransportMatrix, FNoiseMatrix);
    rep->getBackwardJacobianAndNoise(BTransportMatrix, BNoiseMatrix);

    ReferenceStateOnPlane* refState = new ReferenceStateOnPlane(*seedState);
    fitterInfo->setReferenceState(refState);

    if (i==0) { // if we are at first measurement and seed state is defined somewhere else, still set forward info to default
      segmentLen = 0;
    }

    refState->setForwardSegmentLength(segmentLen);
    refState->setForwardTransportMatrix(FTransportMatrix);
    refState->setForwardNoiseMatrix(FNoiseMatrix);

    if (prevFitterInfo != NULL) {
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

    // get MeasurementsOnPlane
    std::vector<AbsMeasurement*> rawMeasurements = trackPoint->getRawMeasurements();
    for ( std::vector< genfit::AbsMeasurement* >::iterator measurement = rawMeasurements.begin(), lastMeasurement =rawMeasurements.end(); measurement != lastMeasurement; ++measurement)
    {
     assert((*measurement) != NULL);
     fitterInfo->setMeasurementsOnPlane((*measurement)->constructMeasurementsOnPlane(rep, plane));
     if (KalmanFiAvailable) {
       fitterInfo->setWeights(oldWeights);
     }
    }

  }

  // set backward info for last reference state
  ReferenceStateOnPlane* prevRefState =  prevFitterInfo->getReferenceState();
  prevRefState->setBackwardSegmentLength(0);
  BTransportMatrix.UnitMatrix();
  prevRefState->setBackwardTransportMatrix(BTransportMatrix);
  BNoiseMatrix.Zero();
  prevRefState->setBackwardNoiseMatrix(BNoiseMatrix);

  // self check
  assert(tr->checkConsistency());
  assert(isTrackPrepared(tr, rep));
}


void
KalmanFitterRefTrack::processTrackPoint(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi, double& chi2, double& ndf, int direction)
{

  unsigned int dim = fi->getRep()->getDim();

  TVectorD dp(dim); // \delta p_{k|k-1}
  TMatrixDSym C(dim); // C_{k|k-1}

  // predict
  if (prevFi != NULL) {
    const TMatrixD& F = fi->getReferenceState()->getTransportMatrix(direction); // Transport matrix
    const TMatrixDSym& N = fi->getReferenceState()->getNoiseMatrix(direction); // Noise matrix
    dp = F * (prevFi->getUpdate(direction)->getState() - prevFi->getReferenceState()->getState());
    C = prevFi->getUpdate(direction)->getCov();
    C.Similarity(F);
    C += N;
    fi->setPrediction(new MeasuredStateOnPlane(dp + fi->getReferenceState()->getState(), C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep()), direction);
#ifdef DEBUG
    std::cout << "\033[31m";
    std::cout << "F (Transport Matrix) "; F.Print();
    std::cout << "Δp_{k-1,k-1} "; (prevFi->getUpdate(direction)->getState() - prevFi->getReferenceState()->getState()).Print();
    std::cout << " p_{k-1,r} (reference state from previous hit)"; prevFi->getReferenceState()->getState().Print();
#endif
  }
  else {
    dp = fi->getPrediction(direction)->getState() - fi->getReferenceState()->getState();
    C = fi->getPrediction(direction)->getCov();
#ifdef DEBUG
    std::cout << "\033[31m";
    std::cout << "p_{k,r} (reference state)"; fi->getReferenceState()->getState().Print();
#endif
  }

#ifdef DEBUG
  std::cout << "Δp_{k|k-1} "; dp.Print();
  std::cout << " p_{k|k-1} (state prediction)"; fi->getPrediction(direction)->getState().Print();
  std::cout << " C_{k|k-1} (covariance prediction)"; C.Print();
  std::cout << "\033[0m";
#endif

  // update
  const MeasurementOnPlane& m = getMeasurement(fi, direction);
  const TVectorD& dm = m.getState() - (m.getHMatrix() * fi->getReferenceState()->getState()); // \delta m_k = m_k - H_k p_{k,r}

  TMatrixDSym covSumInv(C); // (V_k + H_k C_{k|k-1} H_k^T)^(-1)
  covSumInv.Similarity(m.getHMatrix());
  covSumInv += m.getCov();
  tools::invertMatrix(covSumInv);

  TMatrixD CHt(C, TMatrixD::kMultTranspose, m.getHMatrix());

  const TVectorD& res = dm - m.getHMatrix()*dp; //
#ifdef DEBUG
  std::cout << "\033[34m";
  std::cout << "m (measurement) "; m.getState().Print();
  std::cout << "residual        "; res.Print();
  std::cout << "\033[0m";
#endif
  TVectorD updated(TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
  updated += dp;
  updated += fi->getReferenceState()->getState();

  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  C -= covSumInv; // updated Cov

#ifdef DEBUG
  std::cout << " C update "; covSumInv.Print();
  std::cout << "\033[32m";
  std::cout << "Δp_{k|k} "; (updated - fi->getReferenceState()->getState()).Print();
  std::cout << " p_{k|k} (updated state)"; updated.Print();
  std::cout << " C_{k|k} (updated covariance)"; C.Print();
  std::cout << "\033[0m";
#endif

  // Calculate chi² increment.  At the first point chi2inc == 0 and
  // the matrix will not be invertible.
  double chi2inc = 0;
  if (ndf != 0) {
    TMatrixDSym Rinv(C);
    Rinv.Similarity(m.getHMatrix());
    Rinv -= m.getCov();
    Rinv *= -1;
    tools::invertMatrix(Rinv);

    TVectorD resNew(m.getState() - m.getHMatrix()*updated);

    chi2inc = Rinv.Similarity(resNew);
  }
  chi2 += chi2inc;

  double ndfInc = m.getState().GetNrows() * m.getWeight();
  ndf += ndfInc;

#ifdef DEBUG
  std::cout << " chi² inc " << chi2inc << "\n";
  std::cout << " ndf inc  " << ndfInc << "\n";
#endif


  KalmanFittedStateOnPlane* upState = new KalmanFittedStateOnPlane(updated, C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), chi2inc, ndfInc);
  upState->setAuxInfo(fi->getReferenceState()->getAuxInfo());
  fi->setUpdate(upState, direction);

  // check
  assert(fi->checkConsistency());

}
