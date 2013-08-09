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
//#define DEBUGMATH


using namespace genfit;


void KalmanFitterRefTrack::fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{

  if (!isTrackPrepared(tr, rep)) {
    Exception exc("KalmanFitterRefTrack::fitTrack ==> track is not properly prepared.",__LINE__,__FILE__);
    throw exc;
  }

  chi2 = 0;
  ndf = -1. * rep->getDim();
  KalmanFitterInfo* prevFi(NULL);

#ifdef DEBUG
  std::cout << tr->getNumPoints() << " TrackPoints with measurements in this track." << std::endl;
#endif

  bool alreadyFitted(!refitAll_);

  for (size_t i = 0; i < tr->getNumPointsWithMeasurement(); ++i) {
      TrackPoint *tp = 0;
      assert(direction == +1 || direction == -1);
      if (direction == +1)
        tp = tr->getPointWithMeasurement(i);
      else if (direction == -1)
        tp = tr->getPointWithMeasurement(-i-1);

      KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));

      if (alreadyFitted && fi->hasUpdate(direction)) {
        #ifdef DEBUG
        std::cout << "TrackPoint " << i << " is already fitted. \n";
        #endif
        prevFi = fi;
        chi2 += fi->getUpdate(direction)->getChiSquareIncrement();
        ndf += fi->getUpdate(direction)->getNdf();
        continue;
      }

      alreadyFitted = false;

      #ifdef DEBUG
      std::cout << " process TrackPoint nr. " << i << "\n";
      #endif
      processTrackPoint(fi, prevFi, chi2, ndf, direction);

      prevFi = fi;
  }
}


void KalmanFitterRefTrack::processTrack(Track* tr, const AbsTrackRep* rep, bool resortHits)
{

  if (tr->getFitStatus(rep) != NULL && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("KalmanFitterRefTrack::processTrack: Cannot process pruned track!", __LINE__,__FILE__);
    throw exc;
  }

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
      std::cout << " KalmanFitterRefTrack::processTrack with rep " << rep << " (id == " << tr->getIdForRep(rep) << ")"<< ", iteration nr. " << nIt << "\n";
      #endif

      // prepare
      if (!prepareTrack(tr, rep, resortHits) && !refitAll_) {
        #ifdef DEBUG
        std::cout << "KalmanFitterRefTrack::processTrack. Track preparation did not change anything!\n";
        #endif
        status->setIsFitted();
        status->setIsFitConverged();
        status->setHasTrackChanged(false);
        status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
        status->setNumIterations(nIt-1);
        status->setForwardChiSqu(chi2FW);
        status->setBackwardChiSqu(chi2BW);
        status->setForwardNdf(std::max(0., ndfFW));
        status->setBackwardNdf(std::max(0., ndfBW));
        return;
      }

      #ifdef DEBUG
      std::cout << "KalmanFitterRefTrack::processTrack. Prepared Track:";
      tr->Print("C");
      //tr->Print();
      #endif

      // resort
      if (resortHits) {
        if (tr->sort()) {
          #ifdef DEBUG
          std::cout << "KalmanFitterRefTrack::processTrack. Resorted Track:"; tr->Print("C");
          #endif
          prepareTrack(tr, rep, resortHits);// re-prepare if order of hits has changed!
          #ifdef DEBUG
          std::cout << "KalmanFitterRefTrack::processTrack. Prepared resorted Track:"; tr->Print("C");
          #endif
        }
      }


      // fit forward
      #ifdef DEBUG
      std::cout << "KalmanFitterRefTrack::forward fit\n";
      #endif
      fitTrack(tr, rep, chi2FW, ndfFW, +1);

      // fit backward
      #ifdef DEBUG
      std::cout << "KalmanFitterRefTrack::backward fit\n";
      #endif

      // backward fit must not necessarily start at last hit, but if it does, set prediction = forward update and blow up cov
      KalmanFitterInfo* lastInfo = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(-1)->getFitterInfo(rep));
      if (! lastInfo->hasBackwardPrediction()) {
        lastInfo->setBackwardPrediction(new MeasuredStateOnPlane(*(lastInfo->getForwardUpdate())));
        lastInfo->getBackwardPrediction()->blowUpCov(blowUpFactor_);  // blow up cov
        #ifdef DEBUG
        std::cout << "blow up cov for backward fit\n";
        #endif
      }

      fitTrack(tr, rep, chi2BW, ndfBW, -1);

      ++nIt;


      double PvalBW = ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW);
      #ifdef DEBUG
      double PvalFW = ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW);
      #endif


      #ifdef DEBUG
      std::cout << "KalmanFitterRefTrack::Track after fit:"; tr->Print("C");

      std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
          << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
      std::cout << "old pVals: " << oldPvalBW << ", " << oldPvalFW
          << " new pVals: " << PvalBW << ", " << PvalFW << std::endl;
      #endif

      // See if p-value only changed little.  If the initial
      // parameters are very far off, initial chi^2 and the chi^2
      // after the first iteration will be both very close to zero, so
      // we need to have at least two iterations here.  Convergence
      // doesn't make much sense before running twice anyway.
      if (nIt > 1 && fabs(oldPvalBW - PvalBW) < deltaPval_)  {
        // Finished
        #ifdef DEBUG
        std::cout << "Fit is converged! \n";
        #endif
        status->setIsFitConverged();
        break;
      }
      else {
        oldPvalBW = PvalBW;
        #ifdef DEBUG
        oldPvalFW = PvalFW;
        oldChi2BW = chi2BW;
        oldChi2FW = chi2FW;
        #endif
      }

      if (nIt >= maxIterations_) {
        #ifdef DEBUG
        std::cout << "KalmanFitterRefTrack::number of max iterations reached!\n";
        #endif
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


  status->setIsFitted();
  status->setHasTrackChanged(false);
  status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
  status->setNumIterations(nIt);
  status->setForwardChiSqu(chi2FW);
  status->setBackwardChiSqu(chi2BW);
  status->setForwardNdf(ndfFW);
  status->setBackwardNdf(ndfBW);

}


bool KalmanFitterRefTrack::prepareTrack(Track* tr, const AbsTrackRep* rep, bool setSortingParams) {

#ifdef DEBUG
  std::cout << "KalmanFitterRefTrack::prepareTrack \n";
#endif

  int notChangedUntil, notChangedFrom;

  // remove outdated reference states
  bool changedSmthg = removeOutdated(tr, rep,  notChangedUntil, notChangedFrom);


  // declare matrices
  TMatrixD FTransportMatrix(rep->getDim(), rep->getDim());
  FTransportMatrix.UnitMatrix();
  TMatrixD BTransportMatrix(rep->getDim(), rep->getDim());

  TMatrixDSym FNoiseMatrix(rep->getDim());
  TMatrixDSym BNoiseMatrix(rep->getDim());

  TVectorD forwardDeltaState(rep->getDim());
  TVectorD backwardDeltaState(rep->getDim());

  // declare stuff
  KalmanFitterInfo* prevFitterInfo(NULL);
  MeasuredStateOnPlane* firstBackwarUpdate(NULL);

  ReferenceStateOnPlane* referenceState(NULL);
  ReferenceStateOnPlane* prevReferenceState(NULL);

  const MeasuredStateOnPlane* smoothedState(NULL);
  const MeasuredStateOnPlane* prevSmoothedState(NULL);

  double trackLen(0);

  bool newRefState(false);

  unsigned int nPoints = tr->getNumPoints();


  unsigned int i=0;

  try {

    // loop over TrackPoints
    for (; i<nPoints; ++i) {

      TrackPoint* trackPoint = tr->getPoint(i);

      // check if we have a measurement
      if (!trackPoint->hasRawMeasurements()) {
        #ifdef DEBUG
        std::cout << "TrackPoint has no rawMeasurements -> continue \n";
        #endif
        continue;
      }


      // get fitterInfo
      KalmanFitterInfo* fitterInfo(NULL);
      if (trackPoint->hasFitterInfo(rep))
        fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

      // create new fitter info if none available
      if (fitterInfo == NULL) {
        #ifdef DEBUG
        std::cout << "create new KalmanFitterInfo \n";
        #endif
        fitterInfo = new KalmanFitterInfo(trackPoint, rep);
        trackPoint->setFitterInfo(fitterInfo);
        changedSmthg = true;
      }
      else {
        #ifdef DEBUG
        std::cout << "TrackPoint " << i << " (" << trackPoint << ") already has KalmanFitterInfo \n";
        #endif

        if (prevFitterInfo == NULL) {
          if (fitterInfo->hasBackwardUpdate()) {
            firstBackwarUpdate = new MeasuredStateOnPlane(*(fitterInfo->getBackwardUpdate()));
            SOPsToDestruct_.push_back(firstBackwarUpdate);
          }

        }
      }

      // get smoothedState if available
      if (fitterInfo->hasPredictionsAndUpdates()) {
        smoothedState = fitterInfo->getFittedState(true);
      }
      else {
        smoothedState = NULL;
      }


      // if fitterInfo already has a reference state, continue
      if (fitterInfo->hasReferenceState()) {

        referenceState = fitterInfo->getReferenceState();
        prevFitterInfo = fitterInfo;
        prevSmoothedState = smoothedState;

        if (!newRefState) {
          #ifdef DEBUG
          std::cout << "TrackPoint already has referenceState and previous referenceState has not been altered -> continue \n";
          #endif
          prevReferenceState = referenceState;
          trackLen += referenceState->getForwardSegmentLength();
          if (setSortingParams)
            trackPoint->setSortingParameter(trackLen);
          continue;
        }

        // previous refState has been altered ->need to update transport matrices
        #ifdef DEBUG
        std::cout << "TrackPoint already has referenceState but previous referenceState has been altered -> update transport matrices and continue \n";
        #endif
        StateOnPlane* stateToExtrapolate = new StateOnPlane(*prevReferenceState);
        SOPsToDestruct_.push_back(stateToExtrapolate);

		// make sure track is consistent if extrapolation fails
        prevReferenceState->resetBackward();
        referenceState->resetForward();

        double segmentLen = rep->extrapolateToPlane(*stateToExtrapolate, fitterInfo->getReferenceState()->getPlane(), false, true);
        #ifdef DEBUG
        std::cout << "extrapolated stateToExtrapolate (prevReferenceState) by " << segmentLen << " cm.\n";
        #endif
        trackLen += segmentLen;

        rep->getForwardJacobianAndNoise(FTransportMatrix, FNoiseMatrix, forwardDeltaState);
        rep->getBackwardJacobianAndNoise(BTransportMatrix, BNoiseMatrix, backwardDeltaState);

        prevReferenceState->setBackwardSegmentLength(-segmentLen);
        prevReferenceState->setBackwardTransportMatrix(BTransportMatrix);
        prevReferenceState->setBackwardNoiseMatrix(BNoiseMatrix);
        prevReferenceState->setBackwardDeltaState(backwardDeltaState);

        referenceState->setForwardSegmentLength(segmentLen);
        referenceState->setForwardTransportMatrix(FTransportMatrix);
        referenceState->setForwardNoiseMatrix(FNoiseMatrix);
        referenceState->setForwardDeltaState(forwardDeltaState);

        if (setSortingParams)
          trackPoint->setSortingParameter(trackLen);


        prevReferenceState = referenceState;
        newRefState = false;

        continue;
      }

      newRefState = false;


      // Construct plane
      SharedPlanePtr plane;
      if (smoothedState != NULL) {
        #ifdef DEBUG
        std::cout << "construct plane with smoothedState \n";
        #endif
        plane = trackPoint->getRawMeasurement(0)->constructPlane(smoothedState);
      }
      else if (prevSmoothedState != NULL) {
        #ifdef DEBUG
        std::cout << "construct plane with prevSmoothedState \n";
        #endif
        plane = trackPoint->getRawMeasurement(0)->constructPlane(prevSmoothedState);
      }
      else if (prevReferenceState != NULL) {
        #ifdef DEBUG
        std::cout << "construct plane with prevReferenceState \n";
        #endif
        plane = trackPoint->getRawMeasurement(0)->constructPlane(prevReferenceState);
      }
      else if (rep != tr->getCardinalRep() &&
                dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
                static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
        #ifdef DEBUG
        std::cout << "construct plane with smoothed state of cardinal rep fit \n";
        #endif
        TVector3 pos, mom;
        tr->getCardinalRep()->getPosMom(*static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true), pos, mom);
        StateOnPlane* cardinalState = new StateOnPlane(rep);
        SOPsToDestruct_.push_back(cardinalState);
        rep->setPosMom(*cardinalState, pos, mom); // also fills auxInfo
        plane = trackPoint->getRawMeasurement(0)->constructPlane(cardinalState);
      }
      else {
        #ifdef DEBUG
        std::cout << "construct plane with state from track \n";
        #endif
        StateOnPlane* seedFromTrack = new StateOnPlane(rep);
        SOPsToDestruct_.push_back(seedFromTrack);
        rep->setPosMom(*seedFromTrack, tr->getStateSeed()); // also fills auxInfo
        plane = trackPoint->getRawMeasurement(0)->constructPlane(seedFromTrack);
      }

      assert (plane.get() != NULL);



      // do extrapolation and set reference state infos
      StateOnPlane* stateToExtrapolate(NULL);
      if (prevFitterInfo == NULL) { // first measurement
        #ifdef DEBUG
        std::cout << "prevFitterInfo == NULL \n";
        #endif
        if (smoothedState != NULL) {
          #ifdef DEBUG
          std::cout << "extrapolate smoothedState to plane\n";
          #endif
          stateToExtrapolate = new StateOnPlane(*smoothedState);
        }
        else if (referenceState != NULL) {
          #ifdef DEBUG
          std::cout << "extrapolate referenceState to plane\n";
          #endif
          stateToExtrapolate = new StateOnPlane(*referenceState);
        }
        else if (rep != tr->getCardinalRep() &&
                  dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
                  static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
          #ifdef DEBUG
          std::cout << "extrapolate smoothed state of cardinal rep fit to plane\n";
          #endif
          TVector3 pos, mom;
          tr->getCardinalRep()->getPosMom(*static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true), pos, mom);
          stateToExtrapolate = new StateOnPlane(rep);
          rep->setPosMom(*stateToExtrapolate, pos, mom);
        }
        else {
          #ifdef DEBUG
          std::cout << "extrapolate seed from track to plane\n";
          #endif
          stateToExtrapolate = new StateOnPlane(rep);
          rep->setPosMom(*stateToExtrapolate, tr->getStateSeed());
        }
      }
      else {
        assert (prevReferenceState != NULL);
        #ifdef DEBUG
        std::cout << "extrapolate prevReferenceState to plane\n";
        #endif
        stateToExtrapolate = new StateOnPlane(*prevReferenceState);
      }
      SOPsToDestruct_.push_back(stateToExtrapolate);

      // make sure track is consistent if extrapolation fails
      if (prevReferenceState != NULL)
        prevReferenceState->resetBackward();
      fitterInfo->deleteReferenceInfo();

      double segmentLen = rep->extrapolateToPlane(*stateToExtrapolate, plane, false, true);
      trackLen += segmentLen;
      #ifdef DEBUG
      std::cout << "extrapolated stateToExtrapolate by " << segmentLen << " cm.\n";
      #endif

      if (i==0) {
        // if we are at first measurement and seed state is defined somewhere else
        segmentLen = 0;
        trackLen = 0;
      }

      if (setSortingParams)
        trackPoint->setSortingParameter(trackLen);


      // get jacobians and noise matrices
      if (i>0) rep->getForwardJacobianAndNoise(FTransportMatrix, FNoiseMatrix, forwardDeltaState);
      rep->getBackwardJacobianAndNoise(BTransportMatrix, BNoiseMatrix, backwardDeltaState);


      // set backward matrices for previous reference state
      if (prevReferenceState != NULL) {
        prevReferenceState->setBackwardSegmentLength(-segmentLen);
        prevReferenceState->setBackwardTransportMatrix(BTransportMatrix);
        prevReferenceState->setBackwardNoiseMatrix(BNoiseMatrix);
        prevReferenceState->setBackwardDeltaState(backwardDeltaState);
      }


      // create new reference state
      newRefState = true;
      referenceState = new ReferenceStateOnPlane(stateToExtrapolate->getState(),
						 stateToExtrapolate->getPlane(),
						 stateToExtrapolate->getRep(),
						 stateToExtrapolate->getAuxInfo());
      referenceState->setForwardSegmentLength(segmentLen);
      referenceState->setForwardTransportMatrix(FTransportMatrix);
      referenceState->setForwardNoiseMatrix(FNoiseMatrix);
      referenceState->setForwardDeltaState(forwardDeltaState);

      referenceState->resetBackward();

      fitterInfo->setReferenceState(referenceState);


      // get MeasurementsOnPlane
      std::vector<double> oldWeights = fitterInfo->getWeights();
      std::vector<AbsMeasurement*> rawMeasurements = trackPoint->getRawMeasurements();
      for ( std::vector< genfit::AbsMeasurement* >::iterator measurement = rawMeasurements.begin(), lastMeasurement = rawMeasurements.end(); measurement != lastMeasurement; ++measurement) {
        assert((*measurement) != NULL);
        if (measurement == rawMeasurements.begin())
          fitterInfo->setMeasurementsOnPlane((*measurement)->constructMeasurementsOnPlane(rep, plane));
        else
          fitterInfo->addMeasurementsOnPlane((*measurement)->constructMeasurementsOnPlane(rep, plane));
      }
      if (oldWeights.size() == fitterInfo->getNumMeasurements()) {
        fitterInfo->setWeights(oldWeights);
      }

      changedSmthg = true;

      prevReferenceState = referenceState;
      prevFitterInfo = fitterInfo;
      prevSmoothedState = smoothedState;

    } // end loop over TrackPoints

  }
  catch (Exception& e) {
    std::cerr << e.what();

    #ifdef DEBUG
    std::cout << "exception at hit " << i << "\n";
    #endif

    // clean up
    cleanSOPsToDestruct();
    removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

    //prevReferenceState->resetForward();
    //referenceState->resetBackward();

    Exception exc("KalmanFitterRefTrack::prepareTrack: got an exception.",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }




  removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

  if (firstBackwarUpdate != NULL) {
    KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep));
    if (! fi->hasForwardPrediction()) {
      #ifdef DEBUG
      std::cout << "set backwards update of first point as forward prediction (with blown up cov) \n";
      #endif
      if (fi->getPlane() != firstBackwarUpdate->getPlane()) {
        rep->extrapolateToPlane(*firstBackwarUpdate, fi->getPlane());
      }
      firstBackwarUpdate->blowUpCov(blowUpFactor_);
      fi->setForwardPrediction(firstBackwarUpdate);
      SOPsToDestruct_.erase( std::find(SOPsToDestruct_.begin(), SOPsToDestruct_.end(), firstBackwarUpdate) );
    }
  }

  KalmanFitStatus* fitStatus = dynamic_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
  if (fitStatus != NULL)
    fitStatus->setTrackLen(trackLen);

  #ifdef DEBUG
  std::cout << "trackLen of reference track = " << trackLen << "\n";
  #endif


  // clean up
  cleanSOPsToDestruct();


  // self check
  //assert(tr->checkConsistency());
  assert(isTrackPrepared(tr, rep));

  return changedSmthg;
}


bool
KalmanFitterRefTrack::removeOutdated(Track* tr, const AbsTrackRep* rep, int& notChangedUntil, int& notChangedFrom) const {

  #ifdef DEBUG
  std::cout << "KalmanFitterRefTrack::removeOutdated \n";
  #endif

  bool changedSmthg(false);

  unsigned int nPoints = tr->getNumPoints();
  notChangedUntil = nPoints-1;
  notChangedFrom = 0;

  // loop over TrackPoints
  for (unsigned int i=0; i<nPoints; ++i) {

    TrackPoint* trackPoint = tr->getPoint(i);

    // check if we have a measurement
    if (!trackPoint->hasRawMeasurements()) {
      #ifdef DEBUG
      std::cout << "TrackPoint has no rawMeasurements -> continue \n";
      #endif
      continue;
    }

    // get fitterInfo
    KalmanFitterInfo* fitterInfo(NULL);
    if (trackPoint->hasFitterInfo(rep))
      fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

    if (fitterInfo == NULL)
      continue;


    // check if we need to calculate or update reference state
    if (fitterInfo->hasReferenceState()) {

      if (! fitterInfo->hasPredictionsAndUpdates()) {
        #ifdef DEBUG
        std::cout << "reference state but not all predictions & updates -> do not touch reference state. \n";
        #endif
        continue;
      }


      const MeasuredStateOnPlane* smoothedState = fitterInfo->getFittedState(true);
      TVectorD res(smoothedState->getState() - fitterInfo->getReferenceState()->getState());
      double chi2(0);

      // calculate chi2, ignore off diagonals
      for (int j=0; j<smoothedState->getCov().GetNcols(); ++j)
        chi2 += res[j]*res[j] / smoothedState->getCov()(j,j);

      if (chi2 < deltaChi2Ref_) {
        // reference state is near smoothed state ->  do not update reference state
        #ifdef DEBUG
        std::cout << "reference state is near smoothed state ->  do not update reference state, chi2 = " << chi2 << "\n";
        #endif
        continue;
      }
      #ifdef DEBUG
      else {
        std::cout << "reference state is not close to smoothed state, chi2 = " << chi2 << "\n";
      }
      #endif
    }


    #ifdef DEBUG
    std::cout << "remove reference info \n";
    #endif

    fitterInfo->deleteReferenceInfo();
    changedSmthg = true;

    // decided to update reference state -> set notChangedUntil (only once)
    if (notChangedUntil == (int)nPoints-1)
      notChangedUntil = i-1;

    notChangedFrom = i+1;

  } // end loop over TrackPoints


  #ifdef DEBUG
  tr->Print("C");
  #endif

  return changedSmthg;
}


void
KalmanFitterRefTrack::removeForwardBackwardInfo(Track* tr, const AbsTrackRep* rep, int notChangedUntil, int notChangedFrom) const {

  unsigned int nPoints = tr->getNumPoints();

  if (refitAll_) {
    tr->deleteForwardInfo(0, -1, rep);
    tr->deleteBackwardInfo(0, -1, rep);
    return;
  }

  // delete forward/backward info from/to points where reference states have changed
  if (notChangedUntil != (int)nPoints-1) {
    tr->deleteForwardInfo(notChangedUntil+1, -1, rep);
  }
  if (notChangedFrom != 0)
    tr->deleteBackwardInfo(0, notChangedFrom-1, rep);

}


void
KalmanFitterRefTrack::processTrackPoint(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi, double& chi2, double& ndf, int direction)
{

#ifdef DEBUG
  std::cout << " KalmanFitterRefTrack::processTrackPoint " << fi->getTrackPoint() << "\n";
#endif

  unsigned int dim = fi->getRep()->getDim();

  TVectorD p(dim); // p_{k|k-1}
  TMatrixDSym C(dim); // C_{k|k-1}

  // predict
  if (prevFi != NULL) {
    const TMatrixD& F = fi->getReferenceState()->getTransportMatrix(direction); // Transport matrix
    assert(F.GetNcols() == (int)dim);
    const TMatrixDSym& N = fi->getReferenceState()->getNoiseMatrix(direction); // Noise matrix
    p = ( F * prevFi->getUpdate(direction)->getState() ) + fi->getReferenceState()->getDeltaState(direction);
    C = prevFi->getUpdate(direction)->getCov();
    C.Similarity(F);
    C += N;
    fi->setPrediction(new MeasuredStateOnPlane(p, C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo()), direction);
#ifdef DEBUGMATH
    std::cout << "\033[31m";
    std::cout << "F (Transport Matrix) "; F.Print();
    std::cout << "p_{k,r} (reference state) "; fi->getReferenceState()->getState().Print();
    std::cout << "c (delta state) "; fi->getReferenceState()->getDeltaState(direction).Print();
    std::cout << "F*p_{k-1,r} + c "; (F *prevFi->getReferenceState()->getState() + fi->getReferenceState()->getDeltaState(direction)).Print();
#endif
  }
  else {
    if (fi->hasPrediction(direction)) {
      #ifdef DEBUG
      std::cout << "  Use prediction as start \n";
      #endif
      p = fi->getPrediction(direction)->getState();
      C = fi->getPrediction(direction)->getCov();
    }
    else {
      #ifdef DEBUG
      std::cout << "  Use reference state and unit cov as start \n";
      #endif
      p = fi->getReferenceState()->getState();
      C.UnitMatrix();
      fi->setPrediction(new MeasuredStateOnPlane(p, C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo()), direction);
    }
#ifdef DEBUGMATH
    std::cout << "\033[31m";
    std::cout << "p_{k,r} (reference state)"; fi->getReferenceState()->getState().Print();
#endif
  }

#ifdef DEBUGMATH
  std::cout << " p_{k|k-1} (state prediction)"; p.Print();
  std::cout << " C_{k|k-1} (covariance prediction)"; C.Print();
  std::cout << "\033[0m";
#endif

  // update
  const MeasurementOnPlane& m = getMeasurement(fi, direction);

  TMatrixDSym covSumInv(C); // (V_k + H_k C_{k|k-1} H_k^T)^(-1)
  covSumInv.Similarity(m.getHMatrix());
  covSumInv += m.getCov();
  tools::invertMatrix(covSumInv);

  TMatrixD CHt(C, TMatrixD::kMultTranspose, m.getHMatrix());

  const TVectorD& res = m.getState() - m.getHMatrix()*p; //
#ifdef DEBUGMATH
  std::cout << "\033[34m";
  std::cout << "m (measurement) "; m.getState().Print();
  std::cout << "V (measurement covariance) "; m.getCov().Print();
  std::cout << "residual        "; res.Print();
  std::cout << "\033[0m";
#endif
  TVectorD updated(TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
#ifdef DEBUGMATH
  std::cout << "\033[32m";
  std::cout << " update"; updated.Print();
  std::cout << "\033[0m";
#endif
  updated += p;

  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  C -= covSumInv; // updated Cov

#ifdef DEBUGMATH
  //std::cout << " C update "; covSumInv.Print();
  std::cout << "\033[32m";
  std::cout << " p_{k|k} (updated state)"; updated.Print();
  std::cout << " C_{k|k} (updated covariance)"; C.Print();
  std::cout << "\033[0m";
#endif

  // Calculate chi² increment.  At the first point chi2inc == 0 and
  // the matrix will not be invertible.
  double chi2inc = 0;
  TVectorD resNew(m.getState() - m.getHMatrix()*updated);
  #ifdef DEBUGMATH
  std::cout << " resNew "; resNew.Print();
  #endif

  // only calculate chi2inc if res != NULL.
  // If matrix inversion fails, chi2inc = 0
  if (resNew != 0) {
    TMatrixDSym Rinv(C);
    Rinv.Similarity(m.getHMatrix());
    Rinv -= m.getCov();
    Rinv *= -1;

    bool couldInvert(true);
    try {
      tools::invertMatrix(Rinv);
    }
    catch (Exception& e) {
      std::cerr << e.what();
      couldInvert = false;
    }

    if (couldInvert)
      chi2inc = Rinv.Similarity(resNew);
  }

  chi2 += chi2inc;


  double ndfInc = m.getState().GetNrows() * m.getWeight();
  ndf += ndfInc;

#ifdef DEBUG
  std::cout << " chi² inc " << chi2inc << "\n";
  std::cout << " ndf inc  " << ndfInc << "\n";
#endif


  KalmanFittedStateOnPlane* upState = new KalmanFittedStateOnPlane(updated, C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo(), chi2inc, ndfInc);
  upState->setAuxInfo(fi->getReferenceState()->getAuxInfo());
  fi->setUpdate(upState, direction);

  // check
  assert(fi->checkConsistency());

}


void
KalmanFitterRefTrack::cleanSOPsToDestruct() {
  for (std::vector<const StateOnPlane*>::iterator it = SOPsToDestruct_.begin(); it != SOPsToDestruct_.end(); ++it) {
    delete *it;
  }
  SOPsToDestruct_.clear();
}



