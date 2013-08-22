/* Copyright 2013, Ludwig-Maximilians Universität München, Technische Universität München
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

/* This implements the Kalman fitter with reference track.  */

#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"
#include "Exception.h"

#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"

#include "boost/scoped_ptr.hpp"

#include <Math/ProbFunc.h>

namespace {
  const bool debug = false;
  const bool debugMath = false;
}

using namespace genfit;


TrackPoint* KalmanFitterRefTrack::fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{

  //if (!isTrackPrepared(tr, rep)) {
  //  Exception exc("KalmanFitterRefTrack::fitTrack ==> track is not properly prepared.",__LINE__,__FILE__);
  //  throw exc;
  //}

  chi2 = 0;
  ndf = -1. * rep->getDim();
  KalmanFitterInfo* prevFi(NULL);

  TrackPoint* retVal(NULL);

  if (debug)
    std::cout << tr->getNumPoints() << " TrackPoints with measurements in this track." << std::endl;

  bool alreadyFitted(!refitAll_);

  for (size_t i = 0; i < tr->getNumPointsWithMeasurement(); ++i) {
      TrackPoint *tp = 0;
      assert(direction == +1 || direction == -1);
      if (direction == +1)
        tp = tr->getPointWithMeasurement(i);
      else if (direction == -1)
        tp = tr->getPointWithMeasurement(-i-1);

      if (! tp->hasFitterInfo(rep)) {
	if (debug)
	  std::cout << "TrackPoint " << i << " has no fitterInfo -> continue. \n";
        continue;
      }

      KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));

      if (alreadyFitted && fi->hasUpdate(direction)) {
	if (debug)
	  std::cout << "TrackPoint " << i << " is already fitted -> continue. \n";
        prevFi = fi;
        chi2 += fi->getUpdate(direction)->getChiSquareIncrement();
        ndf += fi->getUpdate(direction)->getNdf();
        continue;
      }

      alreadyFitted = false;

      if (debug)
	std::cout << " process TrackPoint nr. " << i << "\n";
      processTrackPoint(fi, prevFi, chi2, ndf, direction);
      retVal = tp;

      prevFi = fi;
  }

  return retVal;
}


void KalmanFitterRefTrack::processTrack(Track* tr, const AbsTrackRep* rep, bool resortHits)
{

  if (tr->getFitStatus(rep) != NULL && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("KalmanFitterRefTrack::processTrack: Cannot process pruned track!", __LINE__,__FILE__);
    throw exc;
  }

  double oldChi2FW = 1e6;
  double oldPvalFW = 0.;
  double oldChi2BW = 1e6;
  double oldPvalBW = 0.;
  double chi2FW(0), ndfFW(0);
  double chi2BW(0), ndfBW(0);

  KalmanFitStatus* status = new KalmanFitStatus();
  tr->setFitStatus(status, rep);

  status->setIsFittedWithReferenceTrack(true);

  unsigned int nIt=0;
  for (;;) {

    try {
      if (debug)
	std::cout << " KalmanFitterRefTrack::processTrack with rep " << rep
		  << " (id == " << tr->getIdForRep(rep) << ")"<< ", iteration nr. " << nIt << "\n";

      // prepare
      if (!prepareTrack(tr, rep, resortHits) && !refitAll_) {
	if (debug)
	  std::cout << "KalmanFitterRefTrack::processTrack. Track preparation did not change anything!\n";
        status->setIsFitted();
        status->setIsFitConverged();
        status->setHasTrackChanged(false);
        status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
        status->setNumIterations(nIt);
        status->setForwardChi2(chi2FW);
        status->setBackwardChi2(chi2BW);
        status->setForwardNdf(std::max(0., ndfFW));
        status->setBackwardNdf(std::max(0., ndfBW));
	if (debug)
	  status->Print();
        return;
      }

      if (debug) {
	std::cout << "KalmanFitterRefTrack::processTrack. Prepared Track:";
	tr->Print("C");
	//tr->Print();
      }

      // resort
      if (resortHits) {
        if (tr->sort()) {
	  if (debug) {
	    std::cout << "KalmanFitterRefTrack::processTrack. Resorted Track:";
	    tr->Print("C");
	  }
          prepareTrack(tr, rep, resortHits);// re-prepare if order of hits has changed!
	  if (debug) {
	    std::cout << "KalmanFitterRefTrack::processTrack. Prepared resorted Track:";
	    tr->Print("C");
	  }
        }
      }


      // fit forward
      if (debug)
	std::cout << "KalmanFitterRefTrack::forward fit\n";
      TrackPoint* lastProcessedPoint = fitTrack(tr, rep, chi2FW, ndfFW, +1);

      // fit backward
      if (debug)
	std::cout << "KalmanFitterRefTrack::backward fit\n";

      // backward fit must not necessarily start at last hit, set prediction = forward update and blow up cov
      if (lastProcessedPoint != NULL) {
        KalmanFitterInfo* lastInfo = static_cast<KalmanFitterInfo*>(lastProcessedPoint->getFitterInfo(rep));
        if (! lastInfo->hasBackwardPrediction()) {
          lastInfo->setBackwardPrediction(new MeasuredStateOnPlane(*(lastInfo->getForwardUpdate())));
          lastInfo->getBackwardPrediction()->blowUpCov(blowUpFactor_);  // blow up cov
	  if (debug)
	    std::cout << "blow up cov for backward fit at TrackPoint " << lastProcessedPoint << "\n";
        }
      }

      fitTrack(tr, rep, chi2BW, ndfBW, -1);

      ++nIt;


      double PvalBW = ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW);
      double PvalFW = debug ? ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW) : 0; // Don't calculate if not debugging as this function potentially takes a lot of time.

      if (debug) {
	std::cout << "KalmanFitterRefTrack::Track after fit:"; tr->Print("C");

	std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
		  << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
	std::cout << "old pVals: " << oldPvalBW << ", " << oldPvalFW
		  << " new pVals: " << PvalBW << ", " << PvalFW << std::endl;
      }

      // See if p-value only changed little.  If the initial
      // parameters are very far off, initial chi^2 and the chi^2
      // after the first iteration will be both very close to zero, so
      // we need to have at least two iterations here.  Convergence
      // doesn't make much sense before running twice anyway.
      bool converged(false);
      bool finished(false);
      if (nIt > 1 && fabs(oldPvalBW - PvalBW) < deltaPval_)  {
        // if pVal ~ 0, check if chi2 has changed a lot
        if (PvalBW < 0.001 &&
            fabs(oldChi2BW / chi2BW) > 2.) {
          finished = false;
        }
        else {
          finished = true;
          converged = true;
        }

        if (PvalBW == 0.)
          converged = false;
      }

      if (finished) {
	if (debug) {
	  std::cout << "Fit is finished! ";
	  if(converged)
	    std::cout << "Fit is converged! ";
	  std::cout << "\n";
	}
        status->setIsFitConverged(converged);
        break;
      }
      else {
        oldPvalBW = PvalBW;
        oldChi2BW = chi2BW;
	if (debug) {
	  oldPvalFW = PvalFW;
	  oldChi2FW = chi2FW;
	}
      }

      if (nIt >= maxIterations_) {
	if (debug)
	  std::cout << "KalmanFitterRefTrack::number of max iterations reached!\n";
        break;
      }
    }
    catch(Exception& e) {
      std::cerr << e.what();
      status->setIsFitted(false);
      status->setIsFitConverged(false);
      if (debug)
	status->Print();
      return;
    }

  }


  TrackPoint* tp = tr->getPointWithMeasurementAndFitterInfo(0, rep);
  if (tp != NULL &&
      static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep))->hasBackwardUpdate())
    status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));

  if (tp != NULL) {
    status->setIsFitted();
  }
  else { // none of the trackPoints has a fitterInfo
    status->setIsFitted(false);
    status->setIsFitConverged(false);
  }

  status->setHasTrackChanged(false);
  status->setNumIterations(nIt);
  status->setForwardChi2(chi2FW);
  status->setBackwardChi2(chi2BW);
  status->setForwardNdf(ndfFW);
  status->setBackwardNdf(ndfBW);

  if (debug)
    status->Print();
}


bool KalmanFitterRefTrack::prepareTrack(Track* tr, const AbsTrackRep* rep, bool setSortingParams) {

  if (debug)
    std::cout << "KalmanFitterRefTrack::prepareTrack \n";

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
  boost::scoped_ptr<MeasuredStateOnPlane> firstBackwardUpdate;

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
	if (debug)
	  std::cout << "TrackPoint has no rawMeasurements -> continue \n";
        continue;
      }


      // get fitterInfo
      KalmanFitterInfo* fitterInfo(NULL);
      if (trackPoint->hasFitterInfo(rep))
        fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

      // create new fitter info if none available
      if (fitterInfo == NULL) {
	if (debug)
	  std::cout << "create new KalmanFitterInfo \n";
        fitterInfo = new KalmanFitterInfo(trackPoint, rep);
        trackPoint->setFitterInfo(fitterInfo);
        changedSmthg = true;
      }
      else {
	if (debug)
	  std::cout << "TrackPoint " << i << " (" << trackPoint << ") already has KalmanFitterInfo \n";

        if (prevFitterInfo == NULL) {
          if (fitterInfo->hasBackwardUpdate())
	    firstBackwardUpdate.reset(new MeasuredStateOnPlane(*(fitterInfo->getBackwardUpdate())));
        }
      }

      // get smoothedState if available
      if (fitterInfo->hasPredictionsAndUpdates()) {
        smoothedState = &(fitterInfo->getFittedState(true));
      }
      else {
        smoothedState = NULL;
      }


      if (fitterInfo->hasReferenceState()) {

        referenceState = fitterInfo->getReferenceState();
        prevFitterInfo = fitterInfo;
        prevSmoothedState = smoothedState;

        if (!newRefState) {
  	  if (debug)
	    std::cout << "TrackPoint already has referenceState and previous referenceState has not been altered -> continue \n";
          prevReferenceState = referenceState;
          trackLen += referenceState->getForwardSegmentLength();
          if (setSortingParams)
            trackPoint->setSortingParameter(trackLen);
          continue;
        }

        // previous refState has been altered ->need to update transport matrices
	if (debug)
	  std::cout << "TrackPoint already has referenceState but previous referenceState has been altered -> update transport matrices and continue \n";
        StateOnPlane stateToExtrapolate(*prevReferenceState);

        // make sure track is consistent if extrapolation fails
        prevReferenceState->resetBackward();
        referenceState->resetForward();

        double segmentLen = rep->extrapolateToPlane(stateToExtrapolate, fitterInfo->getReferenceState()->getPlane(), false, true);
	if (debug)
	  std::cout << "extrapolated stateToExtrapolate (prevReferenceState) by " << segmentLen << " cm.\n";
        trackLen += segmentLen;

        if (segmentLen == 0) {
          FTransportMatrix.UnitMatrix();
          FNoiseMatrix.Zero();
          forwardDeltaState.Zero();
          BTransportMatrix.UnitMatrix();
          BNoiseMatrix.Zero();
          backwardDeltaState.Zero();
        }
        else {
          rep->getForwardJacobianAndNoise(FTransportMatrix, FNoiseMatrix, forwardDeltaState);
          rep->getBackwardJacobianAndNoise(BTransportMatrix, BNoiseMatrix, backwardDeltaState);
        }

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
	if (debug)
	  std::cout << "construct plane with smoothedState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*smoothedState);
      }
      else if (prevSmoothedState != NULL) {
	if (debug)
	  std::cout << "construct plane with prevSmoothedState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*prevSmoothedState);
      }
      else if (prevReferenceState != NULL) {
	if (debug)
	  std::cout << "construct plane with prevReferenceState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*prevReferenceState);
      }
      else if (rep != tr->getCardinalRep() &&
                trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
                dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
                static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
	if (debug)
	  std::cout << "construct plane with smoothed state of cardinal rep fit \n";
        TVector3 pos, mom;
        const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
        tr->getCardinalRep()->getPosMom(fittedState, pos, mom);
        StateOnPlane cardinalState(rep);
        rep->setPosMom(cardinalState, pos, mom);
        rep->setQop(cardinalState, tr->getCardinalRep()->getQop(fittedState));
        plane = trackPoint->getRawMeasurement(0)->constructPlane(cardinalState);
      }
      else {
	if (debug)
	  std::cout << "construct plane with state from track \n";
        StateOnPlane seedFromTrack(rep);
        rep->setPosMom(seedFromTrack, tr->getStateSeed()); // also fills auxInfo
        plane = trackPoint->getRawMeasurement(0)->constructPlane(seedFromTrack);
      }

      assert (plane.get() != NULL);



      // do extrapolation and set reference state infos
      boost::scoped_ptr<StateOnPlane> stateToExtrapolate(NULL);
      if (prevFitterInfo == NULL) { // first measurement
	if (debug)
	  std::cout << "prevFitterInfo == NULL \n";
        if (smoothedState != NULL) {
	  if (debug)
	    std::cout << "extrapolate smoothedState to plane\n";
          stateToExtrapolate.reset(new StateOnPlane(*smoothedState));
        }
        else if (referenceState != NULL) {
	  if (debug)
	    std::cout << "extrapolate referenceState to plane\n";
          stateToExtrapolate.reset(new StateOnPlane(*referenceState));
        }
        else if (rep != tr->getCardinalRep() &&
                  trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
                  dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
                  static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
	  if (debug)
	    std::cout << "extrapolate smoothed state of cardinal rep fit to plane\n";
          TVector3 pos, mom;
          const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
          tr->getCardinalRep()->getPosMom(fittedState, pos, mom);
          stateToExtrapolate.reset(new StateOnPlane(rep));
          rep->setPosMom(*stateToExtrapolate, pos, mom);
          rep->setQop(*stateToExtrapolate, tr->getCardinalRep()->getQop(fittedState));
        }
        else {
	  if (debug)
	    std::cout << "extrapolate seed from track to plane\n";
          stateToExtrapolate.reset(new StateOnPlane(rep));
          rep->setPosMom(*stateToExtrapolate, tr->getStateSeed());
        }
      } // end if (prevFitterInfo == NULL)
      else {
        if (prevSmoothedState != NULL) {
	  if (debug)
	    std::cout << "extrapolate prevSmoothedState to plane \n";
          stateToExtrapolate.reset(new StateOnPlane(*prevSmoothedState));
        }
        else {
          assert (prevReferenceState != NULL);
	  if (debug)
	    std::cout << "extrapolate prevReferenceState to plane \n";
          stateToExtrapolate.reset(new StateOnPlane(*prevReferenceState));
        }
      }

      // make sure track is consistent if extrapolation fails
      if (prevReferenceState != NULL)
        prevReferenceState->resetBackward();
      fitterInfo->deleteReferenceInfo();

      if (prevFitterInfo != NULL) {
        rep->extrapolateToPlane(*stateToExtrapolate, prevFitterInfo->getPlane());
	if (debug)
	  std::cout << "extrapolated stateToExtrapolate to plane of prevFitterInfo (plane could have changed!) \n";
      }

      double segmentLen = rep->extrapolateToPlane(*stateToExtrapolate, plane, false, true);
      trackLen += segmentLen;
      if (debug) {
	std::cout << "extrapolated stateToExtrapolate by " << segmentLen << " cm.\t";
	std::cout << "charge of stateToExtrapolate: " << rep->getCharge(*stateToExtrapolate) << " \n";
      }

      // get jacobians and noise matrices
      if (segmentLen == 0) {
        FTransportMatrix.UnitMatrix();
        FNoiseMatrix.Zero();
        forwardDeltaState.Zero();
        BTransportMatrix.UnitMatrix();
        BNoiseMatrix.Zero();
        backwardDeltaState.Zero();
      }
      else {
        if (i>0)
          rep->getForwardJacobianAndNoise(FTransportMatrix, FNoiseMatrix, forwardDeltaState);
        rep->getBackwardJacobianAndNoise(BTransportMatrix, BNoiseMatrix, backwardDeltaState);
      }


      if (i==0) {
        // if we are at first measurement and seed state is defined somewhere else
        segmentLen = 0;
        trackLen = 0;
      }
      if (setSortingParams)
        trackPoint->setSortingParameter(trackLen);


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

    if (debug) {
      std::cout << "exception at hit " << i << "\n";
      std::cerr << e.what();
    }

    // clean up
    removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

    // set sorting parameters of TrackPoints where no reference state could be calculated
    for (; i<nPoints; ++i) {
      TrackPoint* trackPoint = tr->getPoint(i);

      if (setSortingParams)
        trackPoint->setSortingParameter(trackLen);

      trackPoint->deleteFitterInfo(rep);
    }

    //prevReferenceState->resetForward();
    //referenceState->resetBackward();

    //Exception exc("KalmanFitterRefTrack::prepareTrack: got an exception.",__LINE__,__FILE__);
    //exc.setFatal();
    //throw exc;

    return true;
  }




  removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

  if (firstBackwardUpdate) {
    KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep));
    if (! fi->hasForwardPrediction()) {
      if (debug)
        std::cout << "set backwards update of first point as forward prediction (with blown up cov) \n";
      if (fi->getPlane() != firstBackwardUpdate->getPlane()) {
        rep->extrapolateToPlane(*firstBackwardUpdate, fi->getPlane());
      }
      firstBackwardUpdate->blowUpCov(blowUpFactor_);
      fi->setForwardPrediction(new MeasuredStateOnPlane(*firstBackwardUpdate));
    }
  }

  KalmanFitStatus* fitStatus = dynamic_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
  if (fitStatus != NULL)
    fitStatus->setTrackLen(trackLen);

  if (debug)
    std::cout << "trackLen of reference track = " << trackLen << "\n";

  // self check
  //assert(tr->checkConsistency());
  assert(isTrackPrepared(tr, rep));

  return changedSmthg;
}


bool
KalmanFitterRefTrack::removeOutdated(Track* tr, const AbsTrackRep* rep, int& notChangedUntil, int& notChangedFrom) const {

  if (debug)
    std::cout << "KalmanFitterRefTrack::removeOutdated \n";

  bool changedSmthg(false);

  unsigned int nPoints = tr->getNumPoints();
  notChangedUntil = nPoints-1;
  notChangedFrom = 0;

  // loop over TrackPoints
  for (unsigned int i=0; i<nPoints; ++i) {

    TrackPoint* trackPoint = tr->getPoint(i);

    // check if we have a measurement
    if (!trackPoint->hasRawMeasurements()) {
      if (debug)
	std::cout << "TrackPoint has no rawMeasurements -> continue \n";
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
        if (debug)
	  std::cout << "reference state but not all predictions & updates -> do not touch reference state. \n";
        continue;
      }


      const MeasuredStateOnPlane& smoothedState = fitterInfo->getFittedState(true);
      TVectorD res(smoothedState.getState() - fitterInfo->getReferenceState()->getState());
      double chi2(0);

      // calculate chi2, ignore off diagonals
      for (int j=0; j<smoothedState.getCov().GetNcols(); ++j)
        chi2 += res[j]*res[j] / smoothedState.getCov()(j,j);

      if (chi2 < deltaChi2Ref_) {
        // reference state is near smoothed state ->  do not update reference state
        if (debug)
	  std::cout << "reference state is near smoothed state ->  do not update reference state, chi2 = " << chi2 << "\n";
        continue;
      } else {
        if (debug)
	  std::cout << "reference state is not close to smoothed state, chi2 = " << chi2 << "\n";
      }
    }

    if (debug)
      std::cout << "remove reference info \n";

    fitterInfo->deleteReferenceInfo();
    changedSmthg = true;

    // decided to update reference state -> set notChangedUntil (only once)
    if (notChangedUntil == (int)nPoints-1)
      notChangedUntil = i-1;

    notChangedFrom = i+1;

  } // end loop over TrackPoints


  if (debug)
    tr->Print("C");

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
  if (debug)
    std::cout << " KalmanFitterRefTrack::processTrackPoint " << fi->getTrackPoint() << "\n";

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
    if (debugMath) {
      std::cout << "\033[31m";
      std::cout << "F (Transport Matrix) "; F.Print();
      std::cout << "p_{k,r} (reference state) "; fi->getReferenceState()->getState().Print();
      std::cout << "c (delta state) "; fi->getReferenceState()->getDeltaState(direction).Print();
      std::cout << "F*p_{k-1,r} + c "; (F *prevFi->getReferenceState()->getState() + fi->getReferenceState()->getDeltaState(direction)).Print();
    }
  }
  else {
    if (fi->hasPrediction(direction)) {
      if (debug)
	std::cout << "  Use prediction as start \n";
      p = fi->getPrediction(direction)->getState();
      C = fi->getPrediction(direction)->getCov();
    }
    else {
      if (debug)
	std::cout << "  Use reference state and seed cov as start \n";
      const AbsTrackRep *rep = fi->getReferenceState()->getRep();
      p = fi->getReferenceState()->getState();

      // Convert from 6D covariance of the seed to whatever the trackRep wants.
      TMatrixDSym dummy(p.GetNrows());
      MeasuredStateOnPlane mop(p, dummy, fi->getReferenceState()->getPlane(), rep, fi->getReferenceState()->getAuxInfo());
      TVector3 pos, mom;
      rep->getPosMom(mop, pos, mom);
      rep->setPosMomCov(mop, pos, mom, fi->getTrackPoint()->getTrack()->getCovSeed());
      // Blow up, set.
      mop.blowUpCov(blowUpFactor_);
      fi->setPrediction(new MeasuredStateOnPlane(mop), direction);
      C = mop.getCov();
    }
    if (debugMath) {
      std::cout << "\033[31m";
      std::cout << "p_{k,r} (reference state)"; fi->getReferenceState()->getState().Print();
    }
  }

  if (debugMath) {
    std::cout << " p_{k|k-1} (state prediction)"; p.Print();
    std::cout << " C_{k|k-1} (covariance prediction)"; C.Print();
    std::cout << "\033[0m";
  }

  // update
  const MeasurementOnPlane& m = getMeasurement(fi, direction);

  TMatrixDSym covSumInv(C); // (V_k + H_k C_{k|k-1} H_k^T)^(-1)
  covSumInv.Similarity(m.getHMatrix());
  covSumInv += m.getCov();
  tools::invertMatrix(covSumInv);

  TMatrixD CHt(C, TMatrixD::kMultTranspose, m.getHMatrix());

  const TVectorD& res = m.getState() - m.getHMatrix()*p; //
  if (debugMath) {
    std::cout << "\033[34m";
    std::cout << "m (measurement) "; m.getState().Print();
    std::cout << "V (measurement covariance) "; m.getCov().Print();
    std::cout << "residual        "; res.Print();
    std::cout << "\033[0m";
  }
  TVectorD updated(TMatrixD(CHt, TMatrixD::kMult, covSumInv) * res);
  if (debugMath) {
    std::cout << "\033[32m";
    std::cout << " update"; updated.Print();
    std::cout << "\033[0m";
  }
  updated += p;

  covSumInv.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
  C -= covSumInv; // updated Cov

  if (debugMath) {
    //std::cout << " C update "; covSumInv.Print();
    std::cout << "\033[32m";
    std::cout << " p_{k|k} (updated state)"; updated.Print();
    std::cout << " C_{k|k} (updated covariance)"; C.Print();
    std::cout << "\033[0m";
  }

  // Calculate chi² increment.  At the first point chi2inc == 0 and
  // the matrix will not be invertible.
  double chi2inc = 0;
  TVectorD resNew(m.getState() - m.getHMatrix()*updated);
  if (debugMath) {
    std::cout << " resNew ";
    resNew.Print();
  }

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

    if (couldInvert) {
      if (debugMath) {
        std::cout << " Rinv ";
	Rinv.Print();
      }
      chi2inc = Rinv.Similarity(resNew);
    }
  }

  chi2 += chi2inc;


  double ndfInc = m.getState().GetNrows();
  if (multipleMeasurementHandling_ == weightedAverage) {
    ndfInc *= m.getWeight();
    //std::cout << "weight = " << m.getWeight() << "\n";
  }

  ndf += ndfInc;


  KalmanFittedStateOnPlane* upState = new KalmanFittedStateOnPlane(updated, C, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo(), chi2inc, ndfInc);
  upState->setAuxInfo(fi->getReferenceState()->getAuxInfo());
  fi->setUpdate(upState, direction);


  if (debug) {
    std::cout << " chi² inc " << chi2inc << "\t";
    std::cout << " ndf inc  " << ndfInc << "\t";
    std::cout << " charge of update  " << fi->getRep()->getCharge(*upState) << "\n";
  }

  // check
  assert(fi->checkConsistency());

}
